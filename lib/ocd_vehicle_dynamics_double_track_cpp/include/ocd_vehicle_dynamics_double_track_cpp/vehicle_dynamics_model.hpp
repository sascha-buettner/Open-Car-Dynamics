// Copyright 2026 Simon Sagmeister

#pragma once
#include <math.h>

#include <algorithm>
#include <array>
#include <eigen3/Eigen/Dense>
#include <memory>
#include <string>
#include <tum_types_cpp/common.hpp>
#include <unordered_map>

#include "ocd_aerodynamics_model_base_cpp/concept.hpp"
#include "ocd_tire_model_base_cpp/concept.hpp"
#include "ocd_types_cpp/types.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/logging.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_equations.hpp"
#include "param_management_cpp/base.hpp"
#include "param_management_cpp/param_reference_manager.hpp"
#include "tsl_logger_cpp/reference_logger.hpp"
namespace tam::ocd::vehicle_dynamics
{
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
class VehicleDynamicsDoubleTrackModel
: public tam::ocd::interfaces::VehicleDynamicsModelBase<
    double_track::States::StateEnum::CNT_LENGTH_STATE_VECTOR, double_track::States::StateNames>
{
  using x = tam::ocd::vehicle_dynamics::double_track::States::StateEnum;
  using double_per_wheel_t = tam::types::common::DataPerWheel<double>;

public:
  VehicleDynamicsDoubleTrackModel();
  // Settings inputs
  void set_wheel_speeds(const double_per_wheel_t & wheel_speeds_radps) override;
  void set_steering_rate(double rate) { eqns_.set_steering_rate(rate); }
  void set_steering_angles(const double_per_wheel_t & steering_angle_per_wheel_rad) override;
  void set_external_influences(const types::ExternalInfluences & external_influences) override;
  void set_x_vec(const StateVectorType & x_vec) override;
  void set_road_surface(std::shared_ptr<const types::RoadSurface> surface)
  {
    eqns_.road_surface_ = std::move(surface);
  }
  // Evaluate the model
  void evaluate() override;
  // Getting outputs
  double_per_wheel_t get_wheel_load() const override;
  double_per_wheel_t get_steering_load() const override;
  types::VehicleDynamicsModelOutput get_vehicle_dynamics_output() const override;
  types::TireContactDiagnostics get_tire_contact_diagnostics() const
  {
    types::TireContactDiagnostics out;
    out.tire_spring_force_N = eqns_.imr_.tire_spring_force_N;
    out.dynamic_radius_m = eqns_.imr_.dynamic_tire_radius_m;
    out.node_normal_N = eqns_.imr_.contact_node_normal_N;
    const int z[4] = {x::z_belt_FL_m, x::z_belt_FR_m, x::z_belt_RL_m, x::z_belt_RR_m};
    const int v[4] = {
      x::v_z_belt_FL_mps, x::v_z_belt_FR_mps, x::v_z_belt_RL_mps, x::v_z_belt_RR_mps};
    for (int i = 0; i < 4; ++i) {
      out.belt_displacement_m[i] = eqns_.x_vec_[z[i]];
      out.belt_velocity_mps[i] = eqns_.x_vec_[v[i]];
    }
    return out;
  }
  struct StationaryDiagnostics
  {
    double_per_wheel_t elastic_normal_N{}, compression_m{}, shaft_position_m{}, endstop_force_N{},
      damper_extrapolated{}, third_extension_m{}, arb_extension_m{}, spring_jacobian{},
      damper_jacobian{}, radial_tangent_Npm{}, radial_reference_m{}, conditioned_Fx_N{},
      actual_normal_N{}, camber_out_rad{};
    std::array<double_per_wheel_t, 5> node_normal_N{};
    double aero_Fz_N{}, aero_My_Nm{}, aero_Mx_Nm{};
    double radial_constitutive_residual{}, radial_root_fallback_used{};
  };
  // Non-const signature only because the legacy per-wheel container lacks const indexing.
  // This accessor does not alter the equations, states, parameters or caches.
  StationaryDiagnostics get_stationary_diagnostics()
  {
    StationaryDiagnostics out;
    auto & e = eqns_;
    for (int i = 0; i < 4; ++i) {
      out.elastic_normal_N[i] = -e.imr_.tire_spring_force_N[i];
      out.actual_normal_N[i] = e.imr_.contact_normal_load_N[i];
      out.camber_out_rad[i] = e.p_.camber_out_rad[i] + e.kinematic_camber_rad_[i];
      out.compression_m[i] = e.imr_.suspension_spring_compression_m[i];
      out.shaft_position_m[i] = e.geometry_equivalent_extension_m_[i] /
                                (i < 2 ? e.p_.damper_motion_ratio_f : e.p_.damper_motion_ratio_r);
      out.endstop_force_N[i] = e.imr_.endstop_force_N[i];
      out.damper_extrapolated[i] = e.imr_.damper_curve_speed_extrapolated[i];
      out.third_extension_m[i] = e.vd_output_.suspension.third_equivalent_extension_m[i];
      out.arb_extension_m[i] = e.kinematic_arb_m_[i];
      out.spring_jacobian[i] = e.kinematic_spring_jac_[i];
      out.damper_jacobian[i] = e.geometry_leverage_[i];
      out.radial_tangent_Npm[i] =
        e.p_.radial_lookup_enabled > .5 ? e.radial_k_Npm_[i] : (i < 2 ? e.p_.c_t_f : e.p_.c_t_r);
      if (e.p_.mf62_vertical_enabled > .5) {
        const auto & t = i < 2 ? e.p_.vertical_front : e.p_.vertical_rear;
        const double radius = e.mf62_loaded_radius_m_[i], h = 1e-7;
        const double pressure = e.p_.inflation_absolute_Pa[i] - e.p_.atmosphere_Pa;
        const double gamma = -(e.p_.camber_out_rad[i] + e.kinematic_camber_rad_[i]);
        const double fx =
          e.p_.mf62_contact_force_coupling > .5 ? e.imr_.tire_forces_tire_frame_N[i][0] : 0.;
        const double fy = e.p_.mf62_contact_force_coupling > .5
                            ? e.imr_.tire_forces_tire_frame_N[i][1] * (i % 2 ? -1. : 1.)
                            : 0.;
        out.radial_tangent_Npm[i] =
          (mf62_vertical_force(t, pressure, e.wheel_speeds_radps_[i], gamma, radius - h, fx, fy) -
           mf62_vertical_force(t, pressure, e.wheel_speeds_radps_[i], gamma, radius + h, fx, fy)) /
          (2 * h);
      }
      out.radial_reference_m[i] = e.imr_.tire_spring_initial_compression_m[i];
      out.conditioned_Fx_N[i] = e.p_.radial_condition_fx_N[i];
      for (int node = 0; node < 5; ++node)
        out.node_normal_N[node][i] = e.imr_.contact_node_normal_N[node][i];
    }
    out.aero_Fz_N = e.imr_.aero_force_N.z;
    out.aero_My_Nm = e.imr_.aero_torque_Nm.y;
    out.aero_Mx_Nm = e.imr_.aero_torque_Nm.x;
    out.radial_constitutive_residual = e.radial_iteration_residual_;
    out.radial_root_fallback_used = e.radial_root_fallback_used_;
    return out;
  }
  StateVectorType get_x_vec() const override;
  StateVectorType get_x_dot_vec() const override;

  // Get logger and param manager
  tam::tsl::LoggerAccessInterface::SharedPtr get_logger() const override;
  tam::pmg::MgmtInterface::SharedPtr get_param_manager() const override;

private:
  VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT> eqns_;
  tam::tsl::ReferenceLogger::SharedPtr logger_ = std::make_shared<tam::tsl::ReferenceLogger>();
  pmg::ParamReferenceManager::SharedPtr param_manager_ =
    std::make_shared<pmg::ParamReferenceManager>();
  int64_t param_state_count_prev_ = 0;

  void declare_parameters();
  void register_log_signals();
};
}  // namespace tam::ocd::vehicle_dynamics
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model_impl.hpp"
