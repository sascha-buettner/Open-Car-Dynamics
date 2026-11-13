// Copyright 2026 Simon Sagmeister

#pragma once
#include <math.h>

#include <eigen3/Eigen/Dense>
#include <functional>
#include <map>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>

#include "ocd_drivetrain_model_base_cpp/concept.hpp"
#include "ocd_steering_actuator_model_base_cpp/concept.hpp"
#include "ocd_types_cpp/road_surface.hpp"
#include "ocd_vehicle_dynamics_model_base_cpp/concept.hpp"
#include "ocd_vehicle_model_base_cpp/base_class.hpp"
#include "ocd_vehicle_model_cpp/conditional_auxiliary_input.hpp"
#include "param_management_cpp/base.hpp"
#include "param_management_cpp/param_manager_composer.hpp"
#include "param_management_cpp/param_value_manager.hpp"
#include "tsl_logger_cpp/composer.hpp"
#include "tum_helpers_cpp/numerical.hpp"
namespace tam::ocd
{
template <
  interfaces::concepts::DrivetrainModel DT_MODEL_T,
  interfaces::concepts::SteeringActuatorModel SA_MODEL_T,
  interfaces::concepts::VehicleDynamicsModel VD_MODEL_T>
class VehicleModel : public auxiliary::conditional_compilation::SaHelperVehicleModel<
                       DT_MODEL_T, SA_MODEL_T, typename SA_MODEL_T::AuxiliaryInputType>

{
private:
  double integration_time_s_{0.};
  bool state_input_recording_{false};
  bool replay_boundary_{false};
  double replay_brake_scale_{1.};
  double excitation_trigger_mode_{0.}, recorded_brake_scale_{1.};
  types::ExternalInfluences recorded_external_influences_{};
  tam::tsl::ValueLogger::SharedPtr state_input_logger_ = std::make_shared<tam::tsl::ValueLogger>();
  void record_state_and_inputs();
  double excitation_mode_{0.}, excitation_hz_{7.}, excitation_amplitude_{0.}, excitation_phase_{0.},
    excitation_burst_{0.}, excitation_lever_{1.0};
  double excitation_onset_{-1.}, excitation_release_{-1.};
  bool excitation_braking_{false};

  // region internal models
  VD_MODEL_T vehicle_dynamics_model_;
  // DT_MODEL_T drivetrain_model_;
  // SA_MODEL_T steering_actuator_model_;
  using auxiliary::conditional_compilation::SaHelperVehicleModel<
    DT_MODEL_T, SA_MODEL_T, typename SA_MODEL_T::AuxiliaryInputType>::drivetrain_model_;
  using auxiliary::conditional_compilation::SaHelperVehicleModel<
    DT_MODEL_T, SA_MODEL_T, typename SA_MODEL_T::AuxiliaryInputType>::steering_actuator_model_;
  //

  using VEHICLE_T = VehicleModel<DT_MODEL_T, SA_MODEL_T, VD_MODEL_T>;

  // region state vector management
  static const size_t OFFSET_VD = 0;
  static const size_t OFFSET_DT = VD_MODEL_T::k_state_vector_length;
  static const size_t OFFSET_SA = OFFSET_DT + DT_MODEL_T::k_state_vector_length;

  static const size_t LEN_STATE_VECTOR = VD_MODEL_T::k_state_vector_length +
                                         DT_MODEL_T::k_state_vector_length +
                                         SA_MODEL_T::k_state_vector_length;
  // end region

  // aliasing for simplification
  using state_vector_t = Eigen::Matrix<double, LEN_STATE_VECTOR, 1>;
  using double_per_wheel_t = tam::types::common::DataPerWheel<double>;

  // region input and outputs
  tam::pmg::ParamValueManager::SharedPtr param_manager_;
  tam::pmg::ParamManagerComposer::SharedPtr param_manager_composed_;

  typename DT_MODEL_T::DriverInputType drivetrain_driver_input_;
  typename SA_MODEL_T::DriverInputType steering_actuator_driver_input_;
  // Make feedback optional
  std::conditional_t<
    std::is_void_v<typename DT_MODEL_T::FeedbackType>, std::monostate,
    typename DT_MODEL_T::FeedbackType>
    drivetrain_feedback_;
  std::conditional_t<
    std::is_void_v<typename SA_MODEL_T::FeedbackType>, std::monostate,
    typename SA_MODEL_T::FeedbackType>
    steering_actuator_feedback_;

  double_per_wheel_t wheel_speeds_radps_{0, 0, 0, 0};
  double_per_wheel_t steering_angle_per_wheel_rad_{0, 0, 0, 0};
  double_per_wheel_t drivetrain_load_torque_per_wheel_Nm_{0, 0, 0, 0};
  double_per_wheel_t steering_load_torque_per_wheel_Nm_{0, 0, 0, 0};
  types::VehicleDynamicsModelOutput vehicle_dynamics_output_{};

  types::ExternalInfluences external_influences_{};

  typename types::VehicleModelOutput model_output_{};

  // debugging
  tam::tsl::ValueLogger::SharedPtr debug_container_ = std::make_shared<tam::tsl::ValueLogger>();
  tam::tsl::LoggerComposer::SharedPtr debug_container_composed_ =
    std::make_shared<tam::tsl::LoggerComposer>(
      std::vector<tam::tsl::LoggerAccessInterface::SharedPtr>{debug_container_});
  Eigen::Matrix<double, LEN_STATE_VECTOR, 1> x = state_vector_t::Zero();

  // region functions for structuring
  void declare_parameters();
  // void update_delay_blocks();
  void set_state_vector_of_models(state_vector_t x_);
  state_vector_t ode(double t, state_vector_t x_);

public:
  VehicleModel();
  tam::tsl::LoggerAccessInterface::SharedPtr get_state_input_logger() const
  {
    return state_input_logger_;
  }
  void set_road_surface(std::shared_ptr<const types::RoadSurface> surface)
  {
    vehicle_dynamics_model_.set_road_surface(std::move(surface));
  }
  auto get_tire_contact_diagnostics() const
    requires requires(const VD_MODEL_T & dynamics) { dynamics.get_tire_contact_diagnostics(); }
  {
    return vehicle_dynamics_model_.get_tire_contact_diagnostics();
  }
  // Diagnostic trajectory replay projects selected continuous states at a
  // simulation boundary. It preserves all other states and actuator history.
  // Production simulations do not call this; it is not a force calibration.
  void project_replay_state(const std::map<std::string, double> & values)
  {
    auto candidate = x;
    for (const auto & [name, value] : values) {
      if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite replay state: " + name);
      bool found = false;
      auto apply = [&](const auto & names, const std::string & prefix, std::size_t offset) {
        for (std::size_t i = 0; i < names.size(); ++i) {
          if (name != prefix + std::string(names[i])) continue;
          candidate[offset + i] = value;
          found = true;
          break;
        }
      };
      apply(VD_MODEL_T::StateNamesTrait::value, "vehicle_dynamics.", OFFSET_VD);
      apply(DT_MODEL_T::StateNamesTrait::value, "drivetrain.", OFFSET_DT);
      apply(SA_MODEL_T::StateNamesTrait::value, "steering_actuator.", OFFSET_SA);
      if (!found) throw std::invalid_argument("Unknown replay state: " + name);
    }
    x = candidate;
    set_state_vector_of_models(x);
  }
  // Offline reconstruction only. The complete accepted state and effective
  // inputs are supplied; no integration or prescribed excitation is reapplied.
  void replay_recorded_boundary(const std::map<std::string, double> & record);

  // model step function
  void step() override;
  // drivetrain block
  void set_drivetrain_input(const typename DT_MODEL_T::DriverInputType & input) override
  {
    drivetrain_driver_input_ = input;
  };
  typename DT_MODEL_T::FeedbackType get_drivetrain_feedback() const override
  {
    if constexpr (!std::is_void_v<typename DT_MODEL_T::FeedbackType>) {
      return drivetrain_feedback_;
    }
  };
  // steering actuator block
  void set_steering_input(const typename SA_MODEL_T::DriverInputType & input) override
  {
    steering_actuator_driver_input_ = input;
  };
  typename SA_MODEL_T::FeedbackType get_steering_actuator_feedback() const override
  {
    if constexpr (!std::is_void_v<typename SA_MODEL_T::FeedbackType>) {
      return steering_actuator_feedback_;
    }
  };
  // vehicle dynamics block
  void set_external_influences(const types::ExternalInfluences & input) override
  {
    external_influences_ = input;
  };
  // Conditionally enable auxiliary input functions
  // Cannot be marked with override since it is a template function
  typename types::VehicleModelOutput get_vehicle_model_output() const override
  {
    return model_output_;
  };
  // Debug Output
  tam::tsl::LoggerAccessInterface::SharedPtr get_logger() const override
  {
    return debug_container_composed_;
  };
  // Parameter handling
  tam::pmg::MgmtInterface::SharedPtr get_param_manager() override
  {
    return param_manager_composed_;
  };
  void reset();
};
}  // namespace tam::ocd
#include "vehicle_model_impl.hpp"
