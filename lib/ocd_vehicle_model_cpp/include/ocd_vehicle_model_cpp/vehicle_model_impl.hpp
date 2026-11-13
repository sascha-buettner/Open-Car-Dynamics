// Copyright 2026 Simon Sagmeister
#pragma once
#include <memory>
#include <string>
#include <vector>

#include "ocd_vehicle_model_base_cpp/base_class.hpp"
#include "vehicle_model.hpp"
namespace tam::ocd
{
#define VEHICLE_T VehicleModel<DT_MODEL_T, SA_MODEL_T, VD_MODEL_T>
template <
  interfaces::concepts::DrivetrainModel DT_MODEL_T,
  interfaces::concepts::SteeringActuatorModel SA_MODEL_T,
  interfaces::concepts::VehicleDynamicsModel VD_MODEL_T>
VEHICLE_T::VehicleModel()
{
  param_manager_ = std::make_shared<tam::pmg::ParamValueManager>();
  std::vector<tam::pmg::MgmtInterface::SharedPtr> param_manager_assembly = {
    vehicle_dynamics_model_.get_param_manager(), drivetrain_model_.get_param_manager(),
    steering_actuator_model_.get_param_manager(), param_manager_};
  param_manager_composed_ =
    std::make_shared<tam::pmg::ParamManagerComposer>(param_manager_assembly);
  // Declare parameters
  declare_parameters();
  model_output_.vehicle_dynamics_output.orientation_rad.x = 0.0;
  model_output_.vehicle_dynamics_output.orientation_rad.y = 0.0;
  model_output_.vehicle_dynamics_output.orientation_rad.z = 0.0;
  // <--
  reset();

  // Assemble composed logger
  debug_container_composed_->register_logger(
    vehicle_dynamics_model_.get_logger(), "vehicle_dynamics/");
  debug_container_composed_->register_logger(
    steering_actuator_model_.get_logger(), "steering_actuator/");
  debug_container_composed_->register_logger(drivetrain_model_.get_logger(), "drivetrain/");
}
template <
  interfaces::concepts::DrivetrainModel DT_MODEL_T,
  interfaces::concepts::SteeringActuatorModel SA_MODEL_T,
  interfaces::concepts::VehicleDynamicsModel VD_MODEL_T>
void VEHICLE_T::step()
{
  // Init local variables
  const double h = param_manager_->get_value("integration_step_size_s").as_double();
  const auto substeps = param_manager_->get_value("integration_substeps").as_int();
  if (!std::isfinite(h) || h <= 0. || substeps < 1 || substeps > 8)
    throw std::invalid_argument(
      "Integration requires positive finite dt and 1..8 synchronous substeps");

  state_input_recording_ = param_manager_->get_value("recording.state_input_only").as_bool();
  excitation_trigger_mode_ =
    param_manager_->get_value("diagnostic_excitation.trigger_mode").as_double();
  excitation_mode_ = param_manager_->get_value("diagnostic_excitation.mode").as_double();
  excitation_hz_ = param_manager_->get_value("diagnostic_excitation.frequency_hz").as_double();
  excitation_amplitude_ = param_manager_->get_value("diagnostic_excitation.amplitude").as_double();
  excitation_phase_ = param_manager_->get_value("diagnostic_excitation.phase_rad").as_double();
  excitation_burst_ = param_manager_->get_value("diagnostic_excitation.burst_cycles").as_double();
  excitation_lever_ = param_manager_->get_value("diagnostic_excitation.aero_lever_m").as_double();
  if (excitation_mode_ != 0.) {
    if (!(excitation_hz_ > 0. && excitation_hz_ <= 100. && excitation_amplitude_ >= 0.))
      throw std::invalid_argument("Invalid diagnostic excitation");
    if (excitation_trigger_mode_ == 1.) {
      excitation_onset_ = 0.;
      excitation_release_ = -1.;
    } else if constexpr (requires { drivetrain_driver_input_.brake_pressure_per_wheel_Pa; }) {
      double pressure = 0.;
      for (int i = 0; i < 4; ++i)
        pressure = std::max(pressure, drivetrain_driver_input_.brake_pressure_per_wheel_Pa[i]);
      if (!excitation_braking_ && pressure >= 1.0) {
        excitation_braking_ = true;
        excitation_onset_ = integration_time_s_;
        excitation_release_ = -1.;
      }
      if (excitation_braking_ && pressure < 1.0) {
        excitation_braking_ = false;
        excitation_release_ = integration_time_s_;
      }
    }
  }
  // Log some stuff
  debug_container_->log("integration_step_size_s", h);
  debug_container_->log("integration_substeps", static_cast<double>(substeps));

  // Bind the function for the integrator
  std::function<state_vector_t(double, state_vector_t)> bound_ode =
    std::bind(&VEHICLE_T::ode, this, std::placeholders::_1, std::placeholders::_2);

  // Calculate the new x - Aka do the integration step
  // Every substep advances the complete coupled state and queries contact at
  // each stage. The external control/readout clock remains h (normally 800 Hz).
  const double inner_h = h / static_cast<double>(substeps);
  for (std::int64_t i = 0; i < substeps; ++i)
    x = tam::helpers::numerical::integration_step_DoPri45<state_vector_t>(
      bound_ode, i * inner_h, x, inner_h);

  // Refresh derivatives and outputs at the accepted state. Dormand-Prince's
  // final evaluated stage is not necessarily the accepted state.
  ode(h, x);
  integration_time_s_ += h;

  model_output_.vehicle_dynamics_output = vehicle_dynamics_output_;
  model_output_.wheel_speeds_radps = wheel_speeds_radps_;
  model_output_.drivetrain_load_torque_per_wheel_Nm = drivetrain_load_torque_per_wheel_Nm_;
  model_output_.steering_angle_per_wheel_rad = steering_angle_per_wheel_rad_;
  model_output_.steering_load_torque_per_wheel_Nm = steering_load_torque_per_wheel_Nm_;
  if (state_input_recording_) record_state_and_inputs();
}
template <
  interfaces::concepts::DrivetrainModel DT_MODEL_T,
  interfaces::concepts::SteeringActuatorModel SA_MODEL_T,
  interfaces::concepts::VehicleDynamicsModel VD_MODEL_T>
void VEHICLE_T::set_state_vector_of_models(state_vector_t x_)
{
  // Set the state vectors first for all the models
  // Vehicle dynamics model
  vehicle_dynamics_model_.set_x_vec(x_.segment(OFFSET_VD, VD_MODEL_T::k_state_vector_length));
  drivetrain_model_.set_x_vec(x_.segment(OFFSET_DT, DT_MODEL_T::k_state_vector_length));
  steering_actuator_model_.set_x_vec(x_.segment(OFFSET_SA, SA_MODEL_T::k_state_vector_length));
}
template <
  interfaces::concepts::DrivetrainModel DT_MODEL_T,
  interfaces::concepts::SteeringActuatorModel SA_MODEL_T,
  interfaces::concepts::VehicleDynamicsModel VD_MODEL_T>
VEHICLE_T::state_vector_t VEHICLE_T::ode(double t, state_vector_t x_)
{
  // Silence warning of unused param
  if constexpr (requires { drivetrain_model_.set_time(0.); })
    drivetrain_model_.set_time(integration_time_s_ + t);

  state_vector_t x_dot_ = state_vector_t::Zero();

  set_state_vector_of_models(x_);  // Set state vector of models

  // Evaluate steering actuator
  steering_actuator_model_.set_driver_input(steering_actuator_driver_input_);
  steering_actuator_model_.set_load(steering_load_torque_per_wheel_Nm_);
  steering_actuator_model_.evaluate();
  steering_angle_per_wheel_rad_ = steering_actuator_model_.get_steering_angles();
  if constexpr (!std::is_void_v<typename SA_MODEL_T::FeedbackType>) {
    steering_actuator_feedback_ = steering_actuator_model_.get_feedback();
  }
  x_dot_.segment(OFFSET_SA, SA_MODEL_T::k_state_vector_length) =
    steering_actuator_model_.get_x_dot_vec();

  // Use the current ODE-stage drivetrain state, never the preceding stage's cache.
  wheel_speeds_radps_ = drivetrain_model_.get_wheel_speeds();
  // region eval the vd model
  vehicle_dynamics_model_.set_wheel_speeds(wheel_speeds_radps_);
  vehicle_dynamics_model_.set_steering_angles(steering_angle_per_wheel_rad_);
  if constexpr (requires {
                  steering_actuator_model_.get_mean_steering_rate();
                  vehicle_dynamics_model_.set_steering_rate(0.);
                })
    vehicle_dynamics_model_.set_steering_rate(steering_actuator_model_.get_mean_steering_rate());
  auto applied = external_influences_;
  double added_pitch = 0., added_heave = 0., brake_scale = 1., grip_scale = 1., wave = 0.;
  double_per_wheel_t added_wheel{0, 0, 0, 0};
  if (!replay_boundary_ && excitation_mode_ != 0. && excitation_onset_ >= 0.) {
    const double now = integration_time_s_ + t, age = now - excitation_onset_, ramp = .05;
    double envelope = age <= 0. ? 0. : age < ramp ? .5 * (1. - std::cos(M_PI * age / ramp)) : 1.;
    if (excitation_release_ >= 0.)
      envelope *= std::clamp(1. - (now - excitation_release_) / ramp, 0., 1.);
    if (excitation_burst_ > 0.) {
      const double end = excitation_burst_ / excitation_hz_;
      if (age >= end)
        envelope = 0.;
      else if (age > end - ramp)
        envelope *= .5 * (1. - std::cos(M_PI * (end - age) / ramp));
    }
    const double phase = 2 * M_PI * excitation_hz_ * age + excitation_phase_;
    wave = envelope * std::sin(phase);
    if (excitation_mode_ == 1.) added_pitch = excitation_amplitude_ * wave;
    if (excitation_mode_ == 2.) {
      added_heave = excitation_amplitude_ * wave;
      added_pitch = -excitation_lever_ * added_heave;
    }
    if (excitation_mode_ == 3.) {
      const double speed = std::max(10., std::abs(x_[6]));
      const double length = vehicle_dynamics_model_.get_param_manager()
                              ->get_value("vehicle_dynamics_double_track.wheelbase_m")
                              .as_double();
      added_wheel[0] = added_wheel[1] = .5 * excitation_amplitude_ * wave;
      added_wheel[2] = added_wheel[3] =
        .5 * excitation_amplitude_ * envelope *
        std::sin(phase - 2 * M_PI * excitation_hz_ * length / speed);
      for (int i = 0; i < 4; ++i) applied.tire_force_z_N[i] += added_wheel[i];
    }
    if (excitation_mode_ == 5.) brake_scale = 1. + excitation_amplitude_ * wave;
    if (excitation_mode_ == 6.) {
      grip_scale = 1. + excitation_amplitude_ * wave;
      for (int i = 0; i < 4; ++i) applied.lambda_mue[i] *= grip_scale;
    }
    applied.external_torque_Nm.y += added_pitch;
    applied.external_force_N.z += added_heave;
  }
  if (replay_boundary_) brake_scale = replay_brake_scale_;
  if constexpr (requires { drivetrain_model_.set_external_brake_friction_scale(1.); })
    drivetrain_model_.set_external_brake_friction_scale(brake_scale);
  recorded_external_influences_ = applied;
  recorded_brake_scale_ = brake_scale;
  vehicle_dynamics_model_.set_external_influences(applied);
  debug_container_->log("excitation/native_time_s", integration_time_s_ + t);
  debug_container_->log("excitation/wave", wave);
  debug_container_->log("excitation/pitch_Nm", added_pitch);
  debug_container_->log("excitation/heave_N", added_heave);
  debug_container_->log("excitation/brake_friction_scale", brake_scale);
  debug_container_->log("excitation/tire_friction_scale", grip_scale);
  double power = added_pitch * x_[11] + added_heave * x_[8];
  for (int i = 0; i < 4; ++i) {
    debug_container_->log("excitation/contact_force_N/" + std::to_string(i), added_wheel[i]);
    const auto vd = x_.segment(OFFSET_VD, VD_MODEL_T::k_state_vector_length);
    double velocity = vd[16 + i];
    if (
      vehicle_dynamics_model_.get_param_manager()
        ->get_value("vehicle_dynamics_double_track.tire.contact.belt.enabled")
        .as_double() > .5)
      velocity = vd[36 + i];
    power += added_wheel[i] * velocity;
  }
  debug_container_->log("excitation/mechanical_power_W", power);

  vehicle_dynamics_model_.evaluate();
  drivetrain_load_torque_per_wheel_Nm_ = vehicle_dynamics_model_.get_wheel_load();
  steering_load_torque_per_wheel_Nm_ = vehicle_dynamics_model_.get_steering_load();
  vehicle_dynamics_output_ = vehicle_dynamics_model_.get_vehicle_dynamics_output();
  x_dot_.segment(OFFSET_VD, VD_MODEL_T::k_state_vector_length) =
    vehicle_dynamics_model_.get_x_dot_vec();
  // endregion

  // region eval the drivetrain model
  drivetrain_model_.set_driver_input(drivetrain_driver_input_);
  drivetrain_model_.set_load(drivetrain_load_torque_per_wheel_Nm_);
  drivetrain_model_.evaluate();
  wheel_speeds_radps_ = drivetrain_model_.get_wheel_speeds();
  if constexpr (!std::is_void_v<typename DT_MODEL_T::FeedbackType>) {
    drivetrain_feedback_ = drivetrain_model_.get_feedback();
  }
  x_dot_.segment(OFFSET_DT, DT_MODEL_T::k_state_vector_length) = drivetrain_model_.get_x_dot_vec();
  // endregion

  return x_dot_;
}
template <
  interfaces::concepts::DrivetrainModel DT_MODEL_T,
  interfaces::concepts::SteeringActuatorModel SA_MODEL_T,
  interfaces::concepts::VehicleDynamicsModel VD_MODEL_T>
void VEHICLE_T::reset()
{
  integration_time_s_ = 0.;
  excitation_onset_ = -1.;
  excitation_release_ = -1.;
  excitation_braking_ = false;
  // Reset the state vector with the parameters from the submodels
  for (std::size_t i = 0; i < VD_MODEL_T::k_state_vector_length; i++) {
    x[OFFSET_VD + i] =
      param_manager_
        ->get_value(
          "initial_state.vehicle_dynamics." + std::string(VD_MODEL_T::StateNamesTrait::value[i]))
        .as_double();
  }
  for (std::size_t i = 0; i < DT_MODEL_T::k_state_vector_length; i++) {
    x[OFFSET_DT + i] =
      param_manager_
        ->get_value(
          "initial_state.drivetrain." + std::string(DT_MODEL_T::StateNamesTrait::value[i]))
        .as_double();
  }
  for (std::size_t i = 0; i < SA_MODEL_T::k_state_vector_length; i++) {
    x[OFFSET_SA + i] =
      param_manager_
        ->get_value(
          "initial_state.steering_actuator." + std::string(SA_MODEL_T::StateNamesTrait::value[i]))
        .as_double();
  }

  // Reset the models inputs as well
  steering_actuator_driver_input_ = {};
  drivetrain_driver_input_ = {};
  steering_load_torque_per_wheel_Nm_ = {};
  wheel_speeds_radps_ = {};
  steering_angle_per_wheel_rad_ = {};
  external_influences_ = {};
  drivetrain_load_torque_per_wheel_Nm_ = {};
  vehicle_dynamics_output_ = {};
  model_output_ = {};
  set_state_vector_of_models(x);
}
template <
  interfaces::concepts::DrivetrainModel DT_MODEL_T,
  interfaces::concepts::SteeringActuatorModel SA_MODEL_T,
  interfaces::concepts::VehicleDynamicsModel VD_MODEL_T>
void VEHICLE_T::declare_parameters()
{
  param_manager_->declare_parameter(
    "recording.state_input_only", false, tam::pmg::ParameterType::BOOL,
    "Retain only accepted native state and boundary inputs on the model logging channel");
  param_manager_->declare_parameter(
    "diagnostic_excitation.trigger_mode", 0., tam::pmg::ParameterType::DOUBLE,
    "0=braking onset, 1=whole-run continuous input");
  for (const auto & item : std::vector<std::pair<std::string, double>>{
         {"mode", 0.},
         {"frequency_hz", 7.},
         {"amplitude", 0.},
         {"phase_rad", 0.},
         {"burst_cycles", 0.},
         {"aero_lever_m", .3}})
    param_manager_->declare_parameter(
      "diagnostic_excitation." + item.first, item.second, tam::pmg::ParameterType::DOUBLE,
      "Explicit external experiment; default off");
  // Declare the integration step size for the model
  param_manager_->declare_parameter(
    "integration_step_size_s", 0.0008, tam::pmg::ParameterType::DOUBLE,
    "Integration step size of the model");

  param_manager_->declare_parameter(
    "integration_substeps", std::int64_t{1}, tam::pmg::ParameterType::INTEGER,
    "Synchronous full-state internal integration steps per external simulation tick (1..8)");

  // Declare parameters for the initial states of all submodels
  for (std::size_t i = 0; i < VD_MODEL_T::k_state_vector_length; i++) {
    param_manager_->declare_parameter(
      "initial_state.vehicle_dynamics." + std::string(VD_MODEL_T::StateNamesTrait::value[i]), 0.0,
      tam::pmg::ParameterType::DOUBLE,
      "initial value for the state variable of the vehicle dynamics model");
  }
  for (std::size_t i = 0; i < DT_MODEL_T::k_state_vector_length; i++) {
    param_manager_->declare_parameter(
      std::string("initial_state.drivetrain.") + std::string(DT_MODEL_T::StateNamesTrait::value[i]),
      0.0, tam::pmg::ParameterType::DOUBLE,
      "initial value for the state variable of the drivetrain model");
  }
  for (std::size_t i = 0; i < SA_MODEL_T::k_state_vector_length; i++) {
    param_manager_->declare_parameter(
      "initial_state.steering_actuator." + std::string(SA_MODEL_T::StateNamesTrait::value[i]), 0.0,
      tam::pmg::ParameterType::DOUBLE, "initial value for the state variable of steering actuator");
  }
}
template <
  interfaces::concepts::DrivetrainModel DT_MODEL_T,
  interfaces::concepts::SteeringActuatorModel SA_MODEL_T,
  interfaces::concepts::VehicleDynamicsModel VD_MODEL_T>
void VEHICLE_T::record_state_and_inputs()
{
  auto & logger = *state_input_logger_;
  logger.log("native_time_s", integration_time_s_);
  auto states = [&](const auto & names, const std::string & group, std::size_t offset) {
    for (std::size_t i = 0; i < names.size(); ++i)
      logger.log("state/" + group + "/" + std::string(names[i]), x[offset + i]);
  };
  states(VD_MODEL_T::StateNamesTrait::value, "vehicle_dynamics", OFFSET_VD);
  states(DT_MODEL_T::StateNamesTrait::value, "drivetrain", OFFSET_DT);
  states(SA_MODEL_T::StateNamesTrait::value, "steering_actuator", OFFSET_SA);
  if constexpr (requires { steering_actuator_driver_input_.steering_angle_rad; })
    logger.log("input/steering_angle_rad", steering_actuator_driver_input_.steering_angle_rad);
  if constexpr (requires { drivetrain_driver_input_.throttle; })
    logger.log("input/throttle", drivetrain_driver_input_.throttle);
  if constexpr (requires { drivetrain_driver_input_.push_to_pass; })
    logger.log("input/push_to_pass", double(drivetrain_driver_input_.push_to_pass));
  if constexpr (requires { drivetrain_driver_input_.gear_desired; })
    logger.log("input/gear_desired", drivetrain_driver_input_.gear_desired);
  for (int i = 0; i < 4; ++i) {
    const auto suffix = std::to_string(i);
    if constexpr (requires { drivetrain_driver_input_.brake_pressure_per_wheel_Pa; })
      logger.log(
        "input/brake_pressure_Pa/" + suffix,
        drivetrain_driver_input_.brake_pressure_per_wheel_Pa[i]);
    const std::string wheel_names[4] = {"front_left", "front_right", "rear_left", "rear_right"};
    const auto pressure_key =
      "vehicle_dynamics_double_track.tire.inflation_absolute_Pa." + std::string(wheel_names[i]);
    if (vehicle_dynamics_model_.get_param_manager()->has_parameter(pressure_key))
      logger.log(
        "input/inflation_absolute_Pa/" + suffix,
        vehicle_dynamics_model_.get_param_manager()->get_value(pressure_key).as_double());
    logger.log("input/road_height_m/" + suffix, recorded_external_influences_.z_height_road_m[i]);
    logger.log("input/friction_scale/" + suffix, recorded_external_influences_.lambda_mue[i]);
    logger.log("input/rim_force_z_N/" + suffix, recorded_external_influences_.rim_force_z_N[i]);
    logger.log("input/tire_force_z_N/" + suffix, recorded_external_influences_.tire_force_z_N[i]);
  }
  auto vector_input = [&](const std::string & name, const auto & value) {
    logger.log("input/" + name + "/x", value.x);
    logger.log("input/" + name + "/y", value.y);
    logger.log("input/" + name + "/z", value.z);
  };
  vector_input("external_force_N", recorded_external_influences_.external_force_N);
  vector_input("external_torque_Nm", recorded_external_influences_.external_torque_Nm);
  vector_input("wind_mps", recorded_external_influences_.wind_mps);
  logger.log("input/brake_friction_scale", recorded_brake_scale_);
}
template <
  interfaces::concepts::DrivetrainModel DT_MODEL_T,
  interfaces::concepts::SteeringActuatorModel SA_MODEL_T,
  interfaces::concepts::VehicleDynamicsModel VD_MODEL_T>
void VEHICLE_T::replay_recorded_boundary(const std::map<std::string, double> & record)
{
  auto candidate = x;
  auto require = [&](const std::string & key) {
    const double v = record.at(key);
    if (!std::isfinite(v)) throw std::invalid_argument("Nonfinite recorded boundary: " + key);
    return v;
  };
  // Complete accepted boundaries have a fixed typed layout. Restore each known
  // slot once rather than performing a name search for every state on every tick.
  // The public partial-state projection API and all physical equations stay unchanged.
  auto restore = [&](const auto & names, const std::string & group, std::size_t offset) {
    for (std::size_t i = 0; i < names.size(); ++i)
      candidate[offset + i] = require("state/" + group + "/" + std::string(names[i]));
  };
  restore(VD_MODEL_T::StateNamesTrait::value, "vehicle_dynamics", OFFSET_VD);
  restore(DT_MODEL_T::StateNamesTrait::value, "drivetrain", OFFSET_DT);
  restore(SA_MODEL_T::StateNamesTrait::value, "steering_actuator", OFFSET_SA);
  x = candidate;
  set_state_vector_of_models(x);
  integration_time_s_ = require("native_time_s");
  if constexpr (requires { steering_actuator_driver_input_.steering_angle_rad; })
    steering_actuator_driver_input_.steering_angle_rad = require("input/steering_angle_rad");
  if constexpr (requires { drivetrain_driver_input_.throttle; })
    drivetrain_driver_input_.throttle = require("input/throttle");
  if constexpr (requires { drivetrain_driver_input_.gear_desired; })
    drivetrain_driver_input_.gear_desired = require("input/gear_desired");
  if constexpr (requires { drivetrain_driver_input_.push_to_pass; })
    drivetrain_driver_input_.push_to_pass = require("input/push_to_pass") > .5;
  const std::string wheels[4] = {"front_left", "front_right", "rear_left", "rear_right"};
  for (int i = 0; i < 4; ++i) {
    const auto index = std::to_string(i);
    if constexpr (requires { drivetrain_driver_input_.brake_pressure_per_wheel_Pa; })
      drivetrain_driver_input_.brake_pressure_per_wheel_Pa[i] =
        require("input/brake_pressure_Pa/" + index);
    external_influences_.z_height_road_m[i] = require("input/road_height_m/" + index);
    external_influences_.lambda_mue[i] = require("input/friction_scale/" + index);
    external_influences_.rim_force_z_N[i] = require("input/rim_force_z_N/" + index);
    external_influences_.tire_force_z_N[i] = require("input/tire_force_z_N/" + index);
    const auto key = "vehicle_dynamics_double_track.tire.inflation_absolute_Pa." + wheels[i];
    if (vehicle_dynamics_model_.get_param_manager()->has_parameter(key))
      vehicle_dynamics_model_.get_param_manager()->set_value(
        key, require("input/inflation_absolute_Pa/" + index));
  }
  auto vector = [&](const std::string & name, auto & value) {
    value.x = require("input/" + name + "/x");
    value.y = require("input/" + name + "/y");
    value.z = require("input/" + name + "/z");
  };
  vector("external_force_N", external_influences_.external_force_N);
  vector("external_torque_Nm", external_influences_.external_torque_Nm);
  vector("wind_mps", external_influences_.wind_mps);
  replay_brake_scale_ = require("input/brake_friction_scale");
  if (replay_brake_scale_ < 0.) throw std::domain_error("Negative recorded brake friction scale");
  replay_boundary_ = true;
  try {
    ode(0., x);
  } catch (...) {
    replay_boundary_ = false;
    throw;
  }
  replay_boundary_ = false;
  model_output_.vehicle_dynamics_output = vehicle_dynamics_output_;
  model_output_.wheel_speeds_radps = wheel_speeds_radps_;
  model_output_.drivetrain_load_torque_per_wheel_Nm = drivetrain_load_torque_per_wheel_Nm_;
  model_output_.steering_angle_per_wheel_rad = steering_angle_per_wheel_rad_;
  model_output_.steering_load_torque_per_wheel_Nm = steering_load_torque_per_wheel_Nm_;
}
#undef VEHICLE_T
}  // namespace tam::ocd
