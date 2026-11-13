// Copyright 2026 Simon Sagmeister
#pragma once

#include <string>
#include <vector>

#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
#include "tsl_logger_cpp/type_support.hpp"
namespace tam::ocd::vehicle_dynamics
{
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
tam::tsl::LoggerAccessInterface::SharedPtr
VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::get_logger() const
{
  return logger_;
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
tam::pmg::MgmtInterface::SharedPtr
VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::get_param_manager() const
{
  return param_manager_;
}
// Settings inputs
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::set_wheel_speeds(
  const double_per_wheel_t & wheel_speeds_radps)
{
  eqns_.set_wheel_speeds(wheel_speeds_radps);
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::set_steering_angles(
  const double_per_wheel_t & steering_angle_per_wheel_rad)
{
  eqns_.set_steering_angles(steering_angle_per_wheel_rad);
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::set_external_influences(
  const types::ExternalInfluences & external_influences)
{
  eqns_.set_external_influences(external_influences);
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::set_x_vec(
  const StateVectorType & x_vec)
{
  eqns_.set_x_vec(x_vec);
}
// Getting outputs
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::double_per_wheel_t
VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::get_wheel_load() const
{
  return eqns_.drivetrain_load_torque_per_wheel_Nm_;
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::double_per_wheel_t
VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::get_steering_load() const
{
  return eqns_.steering_load_torque_per_wheel_Nm_;
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
types::VehicleDynamicsModelOutput
VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::get_vehicle_dynamics_output() const
{
  return eqns_.vd_output_;
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>

VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::StateVectorType
VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::get_x_vec() const
{
  return eqns_.x_vec_;
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::StateVectorType
VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::get_x_dot_vec() const
{
  return eqns_.x_dot_vec_;
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::VehicleDynamicsDoubleTrackModel()
{
  declare_parameters();
  register_log_signals();
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::evaluate()
{
  eqns_.evaluate();
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::register_log_signals()
{
  // Assign the debug outputs correctly
  for (int i = 0; i < x::CNT_LENGTH_STATE_VECTOR; i++) {
    logger_->log("x_vec/" + std::string(StateNamesTrait::value[i]), &eqns_.x_vec_[i]);
    logger_->log("x_dot_vec/" + std::string(StateNamesTrait::value[i]), &eqns_.x_dot_vec_[i]);
  }
  logger_->log("imr/road_contact_height_m", &eqns_.imr_.road_contact_height_m);
  for (int node = 0; node < 5; ++node) {
    const auto prefix = "imr/contact_node" + std::to_string(node) + "/";
    logger_->log(prefix + "query_x_m", &eqns_.imr_.contact_node_query_x_m[node]);
    logger_->log(prefix + "query_y_m", &eqns_.imr_.contact_node_query_y_m[node]);
    logger_->log(prefix + "road_z_m", &eqns_.imr_.contact_node_road_z_m[node]);
    logger_->log(prefix + "belt_z_m", &eqns_.imr_.contact_node_belt_z_m[node]);
    logger_->log(prefix + "normal_N", &eqns_.imr_.contact_node_normal_N[node]);
    logger_->log(prefix + "force_x_N", &eqns_.imr_.contact_node_force_x_N[node]);
    logger_->log(prefix + "force_y_N", &eqns_.imr_.contact_node_force_y_N[node]);
  }
  logger_->log("configuration/enveloping_radius_front_m", &eqns_.p_.enveloping_radius_f);
  logger_->log("configuration/enveloping_radius_rear_m", &eqns_.p_.enveloping_radius_r);
  logger_->log("configuration/belt_enabled", &eqns_.p_.belt_enabled);
  logger_->log("configuration/belt_mass_front_kg", &eqns_.p_.belt_mass_f);
  logger_->log("configuration/belt_mass_rear_kg", &eqns_.p_.belt_mass_r);
  logger_->log("imr/belt_carcass_force_N", &eqns_.imr_.belt_carcass_force_N);
  logger_->log("imr/radial_lookup_stiffness_Npm", &eqns_.radial_k_Npm_);
  logger_->log("imr/radial_lookup_preload_N", &eqns_.radial_preload_N_);
  logger_->log("imr/radial_lookup_force_coefficient_Npm", &eqns_.radial_force_coefficient_Npm_);
  logger_->log("imr/radial_lookup_iteration_residual", &eqns_.radial_iteration_residual_);
  logger_->log("imr/radial_root_fallback_used", &eqns_.radial_root_fallback_used_);
  logger_->log("imr/radial_root_function_evaluations", &eqns_.radial_root_evaluations_);

  logger_->log("imr/belt_resulting_force_N", &eqns_.imr_.belt_resulting_force_N);
  logger_->log("imr/belt_elastic_energy_J", &eqns_.imr_.belt_elastic_energy_J);
  logger_->log("imr/rim_support_force_N", &eqns_.imr_.rim_support_force_N);
  logger_->log("configuration/contact_enabled", &eqns_.p_.contact_enabled);
  logger_->log("configuration/strict_surface_coverage", &eqns_.p_.strict_surface_coverage);
  logger_->log("configuration/damper_curve_front_enabled", &eqns_.p_.damper_curve_enabled_f);
  logger_->log("configuration/damper_curve_rear_enabled", &eqns_.p_.damper_curve_enabled_r);
  logger_->log("configuration/heave_stiffness_front_Npm2", &eqns_.p_.heave_stiffness_f);
  logger_->log("configuration/heave_stiffness_rear_Npm2", &eqns_.p_.heave_stiffness_r);
  logger_->log(
    "configuration/height_coupling_deflection_limit_m", &eqns_.p_.lateral_drop_deflection_limit);
  logger_->log("configuration/total_unsprung_front_kg", &eqns_.p_.m_t_f);
  logger_->log("configuration/total_unsprung_rear_kg", &eqns_.p_.m_t_r);
  logger_->log("configuration/sprung_mass_kg", &eqns_.p_dep_.m_sprung);

  logger_->log("configuration/tire_stiffness_front_Npm", &eqns_.p_.c_t_f);
  logger_->log("configuration/tire_stiffness_rear_Npm", &eqns_.p_.c_t_r);
  logger_->log("configuration/cog_height_m", &eqns_.p_.h_cg);
  logger_->log("configuration/lateral_stiffness_front_Npm", &eqns_.p_.lateral_stiffness_f);
  logger_->log("configuration/lateral_stiffness_rear_Npm", &eqns_.p_.lateral_stiffness_r);
  logger_->log("configuration/bump_stiffness_front_Npm2", &eqns_.p_.bump_stiffness_f);
  logger_->log("configuration/bump_stiffness_rear_Npm2", &eqns_.p_.bump_stiffness_r);
  logger_->log("imr/contact_force_x_N", &eqns_.imr_.contact_force_x_N);
  logger_->log("imr/contact_force_y_N", &eqns_.imr_.contact_force_y_N);
  logger_->log("imr/contact_missing_samples", &eqns_.imr_.contact_missing_samples);
  logger_->log("imr/contact_bridged_samples", &eqns_.imr_.contact_bridged_samples);
  logger_->log("imr/tire_lateral_deflection_m", &eqns_.imr_.tire_lateral_deflection_m);
  logger_->log("imr/tire_lateral_drop_m", &eqns_.imr_.tire_lateral_drop_m);
  logger_->log("imr/endstop_force_N", &eqns_.imr_.endstop_force_N);
  logger_->log("imr/suspension_bump_force_N", &eqns_.imr_.suspension_bump_force_N);
  logger_->log("imr/damper_curve_speed_extrapolated", &eqns_.imr_.damper_curve_speed_extrapolated);
  logger_->log("imr/suspension_heave_force_N", &eqns_.imr_.suspension_heave_force_N);
  logger_->log("imr/suspension_heave_energy_J", &eqns_.imr_.suspension_heave_energy_J);
  logger_->log("imr/tire_deformation_extrapolated", &eqns_.imr_.tire_deformation_extrapolated);
  logger_->log("imr/contact_yaw_moment_Nm", &eqns_.imr_.contact_yaw_moment_Nm);
  logger_->log("imr/contact_normal_load_N", &eqns_.imr_.contact_normal_load_N);
  logger_->log(
    "imr/contact_tangential_vertical_force_N", &eqns_.imr_.contact_tangential_vertical_force_N);
  logger_->log("imr/contact_elastic_energy_J", &eqns_.imr_.contact_elastic_energy_J);
  // clang-format off
  logger_->log("external_influences/external_force_N", &eqns_.external_influences_.external_force_N); // NOLINT
  logger_->log("external_influences/external_force_N", &eqns_.external_influences_.external_force_N); // NOLINT
  logger_->log("external_influences/external_torque_Nm", &eqns_.external_influences_.external_torque_Nm); // NOLINT
  logger_->log("external_influences/wind_mps", &eqns_.external_influences_.wind_mps);
  logger_->log("external_influences/rim_force_z_N", &eqns_.external_influences_.rim_force_z_N);
  logger_->log("external_influences/tire_force_z_N", &eqns_.external_influences_.tire_force_z_N);
  logger_->log("external_influences/lambda_mue", &eqns_.external_influences_.lambda_mue);
  logger_->log("external_influences/z_height_road_m", &eqns_.external_influences_.z_height_road_m);
  logger_->log("imr/dynamic_tire_radius", &eqns_.imr_.dynamic_tire_radius_m);
  logger_->log("mf62/free_radius_m",&eqns_.mf62_free_radius_m_);
  logger_->log("mf62/loaded_radius_m",&eqns_.mf62_loaded_radius_m_);
  logger_->log("mf62/deflection_m",&eqns_.mf62_deflection_m_);
  logger_->log("mf62/reference_loaded_radius_m",&eqns_.radial_reference_loaded_m_);
  logger_->log("mf62/gamma_rad",&eqns_.tire_gamma_rad_);
  logger_->log("mf62/normal_iteration_residual_N",&eqns_.radial_residual_N_);
  logger_->log("kinematics/hub_jounce_m",&eqns_.kinematic_q_m_);
  logger_->log("kinematics/spring_equivalent_compression_m",&eqns_.kinematic_spring_m_);
  logger_->log("kinematics/spring_leverage",&eqns_.kinematic_spring_jac_);
  logger_->log("tire/inflation_absolute_Pa",&eqns_.p_.inflation_absolute_Pa);
  logger_->log("kinematics/arb_equivalent_m",&eqns_.kinematic_arb_m_);
  logger_->log("kinematics/arb_leverage",&eqns_.kinematic_arb_jac_);
  logger_->log("tire/overturning_moment_Nm",&eqns_.tire_mx_Nm_);
  logger_->log("tire/self_aligning_moment_Nm",&eqns_.tire_mz_Nm_);
  logger_->log("kinematics/damper_leverage",&eqns_.geometry_leverage_);
  logger_->log("kinematics/contact_x_m",&eqns_.kinematic_contact_x_);
  logger_->log("kinematics/contact_y_m",&eqns_.kinematic_contact_y_);
  logger_->log("kinematics/contact_z_m",&eqns_.kinematic_contact_z_);
  logger_->log("kinematics/longitudinal_jacking_slope",&eqns_.kinematic_jac_x_);
  logger_->log("kinematics/lateral_jacking_slope",&eqns_.kinematic_jac_y_);
  logger_->log("kinematics/rack_m",&eqns_.kinematic_rack_m_);
  logger_->log("kinematics/roll_center_front_lateral_m",&eqns_.roll_center_front_y_m_);
  logger_->log("kinematics/roll_center_rear_lateral_m",&eqns_.roll_center_rear_y_m_);
  logger_->log("geometry/jacking_body_m", &eqns_.kinematic_jacking_m_);
  logger_->log("geometry/camber_delta_rad", &eqns_.kinematic_camber_rad_);
  logger_->log("geometry/toe_delta_rad", &eqns_.kinematic_toe_rad_);
  logger_->log("geometry/roll_center_front_m", &eqns_.p_.h_rc_front);
  logger_->log("geometry/roll_center_rear_m", &eqns_.p_.h_rc_rear);
  logger_->log("imr/effective_steering_angle_per_wheel_rad", &eqns_.imr_.effective_steering_angle_per_wheel_rad); // NOLINT
  logger_->log("imr/tire_spring_initial_compression_m", &eqns_.imr_.tire_spring_initial_compression_m); // NOLINT
  logger_->log("imr/tire_spring_force_N", &eqns_.imr_.tire_spring_force_N);
  logger_->log("imr/vertical_tire_force_N", &eqns_.imr_.vertical_tire_force_N);
  logger_->log("imr/antiroll_bar_force_N", &eqns_.imr_.antiroll_bar_force_N);
  logger_->log("imr/velocity_wheel_over_ground_mps", &eqns_.imr_.velocity_wheel_over_ground_mps);
  logger_->log("imr/velocity_wheel_over_ground_tire_frame_mps", &eqns_.imr_.velocity_wheel_over_ground_tire_frame_mps); // NOLINT
  logger_->log("imr/velocity_tire_rotation_mps", &eqns_.imr_.velocity_tire_rotation_mps);
  logger_->log("imr/tire_longitudinal_slip", &eqns_.imr_.tire_longitudinal_slip);
  logger_->log("imr/tire_slip_angle_rad", &eqns_.imr_.tire_slip_angle_rad);
  logger_->log("imr/tire_forces_tire_frame_N", &eqns_.imr_.tire_forces_tire_frame_N);
  logger_->log("imr/tire_forces_N", &eqns_.imr_.tire_forces_N);
  logger_->log("imr/rolling_resistance_N", &eqns_.imr_.tire_rolling_resistance_N);
  logger_->log("imr/suspension_spring_initial_compression_m", &eqns_.imr_.suspension_spring_initial_compression_m); // NOLINT
  logger_->log("imr/suspension_spring_compression_m", &eqns_.imr_.suspension_spring_compression_m);
  logger_->log("imr/suspension_damper_compression_speed_mps", &eqns_.imr_.suspension_damper_compression_speed_mps); // NOLINT
  logger_->log("imr/suspension_spring_force_N", &eqns_.imr_.suspension_spring_force_N);
  logger_->log("imr/suspension_damper_force_N", &eqns_.imr_.suspension_damper_force_N);
  logger_->log("imr/axle_vertical_force_N", &eqns_.imr_.axle_vertical_force_N);
  logger_->log("imr/resulting_suspension_force_N", &eqns_.imr_.resulting_suspension_force_N);
  logger_->log("imr/resulting_vertical_force_on_wheel_N", &eqns_.imr_.resulting_vertical_force_on_wheel_N); // NOLINT
  logger_->log("imr/aero_force_N", &eqns_.imr_.aero_force_N);
  logger_->log("imr/aero_torque_Nm", &eqns_.imr_.aero_torque_Nm);
  logger_->log("imr/resulting_force_N", &eqns_.imr_.resulting_force_N);
  logger_->log("imr/resulting_torque_Nm", &eqns_.imr_.resulting_torque_Nm);
  // clang-format on

  eqns_.aero_model_.register_log_signals(logger_.get(), "imr/aerodynamics/");

  eqns_.tire_models_.front_left.register_log_signals(logger_.get(), "imr/tire_model/front_left/");
  eqns_.tire_models_.front_right.register_log_signals(logger_.get(), "imr/tire_model/front_right/");
  eqns_.tire_models_.rear_left.register_log_signals(logger_.get(), "imr/tire_model/rear_left/");
  eqns_.tire_models_.rear_right.register_log_signals(logger_.get(), "imr/tire_model/rear_right/");
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>::declare_parameters()
{
  auto p_def_d = [this](
                   std::string param_name, const double & initial_value, double * storage_location,
                   const std::string & description = "") {
    param_manager_->declare_parameter(
      std::string("vehicle_dynamics_double_track.") + param_name, storage_location, initial_value,
      tam::pmg::ParameterType::DOUBLE, description);
  };
  // Masses & Inertial moments
  p_def_d("mass_vehicle_kg", 800, &eqns_.p_.m, "Overall vehicle mass including tires");  // [kg]
  p_def_d("mass_wheel_kg.front", 15, &eqns_.p_.m_t_f, "Mass for an individual tire");    // [kg]
  p_def_d("mass_wheel_kg.rear", 15, &eqns_.p_.m_t_r, "Mass for an individual tire");     // [kg]
  p_def_d("moment_of_inertia_kgpm2.z", 1000, &eqns_.p_.I_z);                             // [kg m^2]
  p_def_d("moment_of_inertia_kgpm2.y", 500, &eqns_.p_.I_y);                              // [kg m^2]
  p_def_d("moment_of_inertia_kgpm2.x", 100, &eqns_.p_.I_x);                              // [kg m^2]

  // Geometric parameters
  p_def_d("wheelbase_m", 3.2, &eqns_.p_.l, "Wheelbase of the vehicle");  // [m]
  p_def_d("track_width_m.front", 1.6, &eqns_.p_.b_f);                    // [m]
  p_def_d("track_width_m.rear", 1.5, &eqns_.p_.b_r);                     // [m]

  // Simulink currently assumes that the cog lies in the axle plane
  p_def_d(
    "cog.distance_from_front_axle_m", 1.724, &eqns_.p_.l_f,
    "Distance from the front axle to the COG");  // [m]
  p_def_d("cog.height_m", 0.3, &eqns_.p_.h_cg);  //  [m]
  p_def_d(
    "suspension.pitch_center_height_m.accel", 0.2, &eqns_.p_.h_pc_accel,
    "Height of the pitch center during acceleration");  //  [m]
  p_def_d(
    "suspension.pitch_center_height_m.decel", 0.2, &eqns_.p_.h_pc_decel,
    "Height of the pitch center during deceleration");  //  [m]
  p_def_d(
    "suspension.roll_center_height_m", 0.2, &eqns_.p_.h_rc, "Height of the roll center");  //  [m]
  p_def_d("suspension.endstop.enabled", 1.0, &eqns_.p_.endstop_enabled);
  p_def_d("suspension.endstop.stiffness_Npm", 1.0, &eqns_.p_.endstop_stiffness);
  p_def_d("suspension.endstop.damping_Nspm", 1.0, &eqns_.p_.endstop_damping);
  p_def_d("suspension.endstop.bump_travel_m.front", 1.0, &eqns_.p_.stroke_bump_front);
  p_def_d("suspension.endstop.rebound_travel_m.front", 1.0, &eqns_.p_.stroke_rebound_front);
  p_def_d("suspension.endstop.bump_travel_m.rear", 1.0, &eqns_.p_.stroke_bump_rear);
  p_def_d("suspension.endstop.rebound_travel_m.rear", 1.0, &eqns_.p_.stroke_rebound_rear);
  p_def_d("suspension.geometry.enabled", 1.0, &eqns_.p_.geometry_enabled);
  p_def_d("suspension.geometry.caster_rad", 1.0, &eqns_.p_.caster_rad);
  p_def_d("suspension.geometry.kingpin_inclination_rad", 1.0, &eqns_.p_.kingpin_inclination_rad);
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.geometry.travel_m", &eqns_.p_.geometry_travel,
    eqns_.p_.geometry_travel, tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Wheel-coordinate geometry lookup; clamped outside measured knots");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.geometry.camber_delta_rad.front",
    &eqns_.p_.geometry_camber_front, eqns_.p_.geometry_camber_front,
    tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Wheel-coordinate geometry lookup; clamped outside measured knots");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.geometry.camber_delta_rad.rear",
    &eqns_.p_.geometry_camber_rear, eqns_.p_.geometry_camber_rear,
    tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Wheel-coordinate geometry lookup; clamped outside measured knots");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.geometry.toe_delta_rad.front",
    &eqns_.p_.geometry_toe_front, eqns_.p_.geometry_toe_front,
    tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Wheel-coordinate geometry lookup; clamped outside measured knots");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.geometry.toe_delta_rad.rear",
    &eqns_.p_.geometry_toe_rear, eqns_.p_.geometry_toe_rear, tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Wheel-coordinate geometry lookup; clamped outside measured knots");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.geometry.roll_center_height_m.front",
    &eqns_.p_.geometry_roll_front, eqns_.p_.geometry_roll_front,
    tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Wheel-coordinate geometry lookup; clamped outside measured knots");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.geometry.roll_center_height_m.rear",
    &eqns_.p_.geometry_roll_rear, eqns_.p_.geometry_roll_rear,
    tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Wheel-coordinate geometry lookup; clamped outside measured knots");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.geometry.steering_mean_angle_rad",
    &eqns_.p_.geometry_steer, eqns_.p_.geometry_steer, tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Wheel-coordinate geometry lookup; clamped outside measured knots");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.geometry.jacking_body_m", &eqns_.p_.geometry_jacking,
    eqns_.p_.geometry_jacking, tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Wheel-coordinate geometry lookup; clamped outside measured knots");

  // Spring & Damping constants
  p_def_d(
    "suspension.vehicle_spring_stiffness_Npm.front", 120000, &eqns_.p_.c_f,
    "Considerung the travel at the virtual suspension point at 0.5 track width");  // [N/m]
  p_def_d(
    "suspension.vehicle_spring_stiffness_Npm.rear", 90000.0, &eqns_.p_.c_r,
    "Considerung the travel at the virtual suspension point at 0.5 track width");  // [N/m]
  p_def_d(
    "suspension.vehicle_damper_coefficient_Nspm.front", 1300.0, &eqns_.p_.d_f,
    "Considerung the travel at the virtual suspension point at 0.5 track width");  // [N*s/m]
  p_def_d(
    "suspension.vehicle_damper_coefficient_Nspm.rear", 900.0, &eqns_.p_.d_r,
    "Considerung the travel at the virtual suspension point at 0.5 track width");  // [N*s/m]
  p_def_d("suspension.damper_curve.enabled.front", 1.0, &eqns_.p_.damper_curve_enabled_f);
  p_def_d("suspension.damper_curve.enabled.rear", 1.0, &eqns_.p_.damper_curve_enabled_r);
  p_def_d(
    "suspension.damper_curve.wheel_damper_motion_ratio.front", 1.0,
    &eqns_.p_.damper_motion_ratio_f);
  p_def_d(
    "suspension.damper_curve.wheel_damper_motion_ratio.rear", 1.0,
    &eqns_.p_.damper_motion_ratio_r);
  for (auto [name, target] : std::initializer_list<std::pair<const char *, std::vector<double> *>>{
         {"speed_mps", &eqns_.p_.damper_speed_mps},
         {"bump_force_N.front", &eqns_.p_.damper_bump_f},
         {"rebound_force_N.front", &eqns_.p_.damper_rebound_f},
         {"bump_force_N.rear", &eqns_.p_.damper_bump_r},
         {"rebound_force_N.rear", &eqns_.p_.damper_rebound_r}})
    param_manager_->declare_parameter(
      std::string("vehicle_dynamics_double_track.suspension.damper_curve.") + name, target, *target,
      tam::pmg::ParameterType::DOUBLE_ARRAY,
      "Measured shaft-speed / shaft-force curve; wheel force uses reciprocal motion ratio");

  p_def_d("tire.atmosphere_Pa", 1.0, &eqns_.p_.atmosphere_Pa);
  p_def_d("tire.radial_lookup.enabled", 1.0, &eqns_.p_.radial_lookup_enabled);
  p_def_d("tire.radial_lookup.robust_solver.enabled", 1.0, &eqns_.p_.radial_robust_solver_enabled);
  p_def_d("tire.low_load_continuation.boundary_N", 1.0, &eqns_.p_.low_load_boundary_N);
  p_def_d("tire.contact.sampling_points", 1.0, &eqns_.p_.contact_sampling_points);
  p_def_d("tire.sidewall.enabled", 1.0, &eqns_.p_.sidewall_enabled);
  p_def_d("tire.mf62_vertical.contact_force_coupling", 1.0, &eqns_.p_.mf62_contact_force_coupling);
  p_def_d("suspension.static_roll_centers.enabled", 1.0, &eqns_.p_.static_roll_centers_enabled);
  p_def_d("suspension.static_roll_centers.height_m.front", 1.0, &eqns_.p_.h_rc_front);
  p_def_d("suspension.static_roll_centers.height_m.rear", 1.0, &eqns_.p_.h_rc_rear);
  p_def_d("tire.load_memory.fraction", 1.0, &eqns_.p_.load_memory_fraction);
  p_def_d("tire.load_memory.build_tau_s", 1.0, &eqns_.p_.load_memory_tau_s);
  p_def_d("tire.fy_radial.mode", 1.0, &eqns_.p_.fy_radial_mode);
  p_def_d("tire.fy_relaxation.strength", 1.0, &eqns_.p_.fy_relaxation_strength);
  p_def_d("tire.force_conditioning.enabled", 1.0, &eqns_.p_.force_conditioning_enabled);
  p_def_d("tire.force_conditioning.kappa.front_left", 1.0, &eqns_.p_.force_condition_kappa[0]);
  p_def_d("tire.force_conditioning.alpha_rad.front_left", 1.0, &eqns_.p_.force_condition_alpha[0]);
  p_def_d("tire.force_conditioning.fy_N.front_left", 1.0, &eqns_.p_.radial_condition_fy_N[0]);
  p_def_d("tire.force_conditioning.kappa.front_right", 1.0, &eqns_.p_.force_condition_kappa[1]);
  p_def_d("tire.force_conditioning.alpha_rad.front_right", 1.0, &eqns_.p_.force_condition_alpha[1]);
  p_def_d("tire.force_conditioning.fy_N.front_right", 1.0, &eqns_.p_.radial_condition_fy_N[1]);
  p_def_d("tire.force_conditioning.kappa.rear_left", 1.0, &eqns_.p_.force_condition_kappa[2]);
  p_def_d("tire.force_conditioning.alpha_rad.rear_left", 1.0, &eqns_.p_.force_condition_alpha[2]);
  p_def_d("tire.force_conditioning.fy_N.rear_left", 1.0, &eqns_.p_.radial_condition_fy_N[2]);
  p_def_d("tire.force_conditioning.kappa.rear_right", 1.0, &eqns_.p_.force_condition_kappa[3]);
  p_def_d("tire.force_conditioning.alpha_rad.rear_right", 1.0, &eqns_.p_.force_condition_alpha[3]);
  p_def_d("tire.force_conditioning.fy_N.rear_right", 1.0, &eqns_.p_.radial_condition_fy_N[3]);

  p_def_d("tire.radial_lookup.stiffness_scale", 1.0, &eqns_.p_.radial_stiffness_scale);
  p_def_d("tire.radial_lookup.nominal_stiffness", 1.0, &eqns_.p_.radial_nominal_stiffness);
  p_def_d("tire.radial_lookup.speed_override.enabled", 1.0, &eqns_.p_.radial_speed_override_enabled);
  p_def_d("tire.radial_lookup.conditioning.enabled", 1.0, &eqns_.p_.radial_conditioning_enabled);
  p_def_d("tire.radial_lookup.conditioning.frozen_preload", 1.0, &eqns_.p_.radial_frozen_preload);
  p_def_d("tire.radial_lookup.conditioning.speed_mps", 1.0, &eqns_.p_.radial_speed_mps);
  p_def_d("tire.radial_lookup.softening_fraction", 1.0, &eqns_.p_.radial_softening_fraction);
  p_def_d("tire.radial_lookup.force_reference_N.front", 1.0, &eqns_.p_.radial_force_reference_f);
  p_def_d("tire.radial_lookup.force_reference_N.rear", 1.0, &eqns_.p_.radial_force_reference_r);
  for (int i = 0; i < 4; ++i) {
    const std::string w =
      std::array<std::string, 4>{"front_left", "front_right", "rear_left", "rear_right"}[i];
    p_def_d(
      "tire.radial_lookup.conditioning.normal_N." + w, 1.0,
      &eqns_.p_.radial_condition_load_N[i]);
    p_def_d(
      "tire.radial_lookup.conditioning.longitudinal_N." + w, 1.0,
      &eqns_.p_.radial_condition_fx_N[i]);
  }

  p_def_d(
    "tire.inflation_absolute_Pa.front_left", 1.0, &eqns_.p_.inflation_absolute_Pa.front_left);
  p_def_d(
    "tire.inflation_absolute_Pa.front_right", 1.0, &eqns_.p_.inflation_absolute_Pa.front_right);
  p_def_d(
    "tire.inflation_absolute_Pa.rear_left", 1.0, &eqns_.p_.inflation_absolute_Pa.rear_left);
  p_def_d(
    "tire.inflation_absolute_Pa.rear_right", 1.0, &eqns_.p_.inflation_absolute_Pa.rear_right);
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.geometry.equivalent_travel_m.front",
    &eqns_.p_.geometry_equiv_front, eqns_.p_.geometry_equiv_front,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Normalized mechanism jounce; derivative1 at nominal");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.geometry.equivalent_travel_m.rear",
    &eqns_.p_.geometry_equiv_rear, eqns_.p_.geometry_equiv_rear,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Normalized mechanism jounce; derivative1 at nominal");
  p_def_d("tire.mf62_vertical.enabled", 1.0, &eqns_.p_.mf62_vertical_enabled);
  p_def_d("tire.mf62_vertical.rolling_enabled", 1.0, &eqns_.p_.mf62_rolling_enabled);
  p_def_d(
    "tire.mf62_vertical.front.UNLOADED_RADIUS", eqns_.p_.vertical_front.UNLOADED_RADIUS,
    &eqns_.p_.vertical_front.UNLOADED_RADIUS);
  p_def_d(
    "tire.mf62_vertical.front.WIDTH", eqns_.p_.vertical_front.WIDTH,
    &eqns_.p_.vertical_front.WIDTH);
  p_def_d(
    "tire.mf62_vertical.front.ASPECT_RATIO", eqns_.p_.vertical_front.ASPECT_RATIO,
    &eqns_.p_.vertical_front.ASPECT_RATIO);
  p_def_d(
    "tire.mf62_vertical.front.RIM_RADIUS", eqns_.p_.vertical_front.RIM_RADIUS,
    &eqns_.p_.vertical_front.RIM_RADIUS);
  p_def_d(
    "tire.mf62_vertical.front.LONGVL", eqns_.p_.vertical_front.LONGVL,
    &eqns_.p_.vertical_front.LONGVL);
  p_def_d(
    "tire.mf62_vertical.front.FNOMIN", eqns_.p_.vertical_front.FNOMIN,
    &eqns_.p_.vertical_front.FNOMIN);
  p_def_d(
    "tire.mf62_vertical.front.NOMPRES", eqns_.p_.vertical_front.NOMPRES,
    &eqns_.p_.vertical_front.NOMPRES);
  p_def_d(
    "tire.mf62_vertical.front.VERTICAL_STIFFNESS", eqns_.p_.vertical_front.VERTICAL_STIFFNESS,
    &eqns_.p_.vertical_front.VERTICAL_STIFFNESS);
  p_def_d(
    "tire.mf62_vertical.front.VERTICAL_DAMPING", eqns_.p_.vertical_front.VERTICAL_DAMPING,
    &eqns_.p_.vertical_front.VERTICAL_DAMPING);
  p_def_d(
    "tire.mf62_vertical.front.Q_RE0", eqns_.p_.vertical_front.Q_RE0,
    &eqns_.p_.vertical_front.Q_RE0);
  p_def_d(
    "tire.mf62_vertical.front.Q_V1", eqns_.p_.vertical_front.Q_V1, &eqns_.p_.vertical_front.Q_V1);
  p_def_d(
    "tire.mf62_vertical.front.Q_V2", eqns_.p_.vertical_front.Q_V2, &eqns_.p_.vertical_front.Q_V2);
  p_def_d(
    "tire.mf62_vertical.front.Q_FZ2", eqns_.p_.vertical_front.Q_FZ2,
    &eqns_.p_.vertical_front.Q_FZ2);
  p_def_d(
    "tire.mf62_vertical.front.Q_FCX", eqns_.p_.vertical_front.Q_FCX,
    &eqns_.p_.vertical_front.Q_FCX);
  p_def_d(
    "tire.mf62_vertical.front.Q_FCY", eqns_.p_.vertical_front.Q_FCY,
    &eqns_.p_.vertical_front.Q_FCY);
  p_def_d(
    "tire.mf62_vertical.front.Q_FCY2", eqns_.p_.vertical_front.Q_FCY2,
    &eqns_.p_.vertical_front.Q_FCY2);
  p_def_d(
    "tire.mf62_vertical.front.Q_CAM1", eqns_.p_.vertical_front.Q_CAM1,
    &eqns_.p_.vertical_front.Q_CAM1);
  p_def_d(
    "tire.mf62_vertical.front.Q_CAM2", eqns_.p_.vertical_front.Q_CAM2,
    &eqns_.p_.vertical_front.Q_CAM2);
  p_def_d(
    "tire.mf62_vertical.front.Q_CAM3", eqns_.p_.vertical_front.Q_CAM3,
    &eqns_.p_.vertical_front.Q_CAM3);
  p_def_d(
    "tire.mf62_vertical.front.Q_FYS1", eqns_.p_.vertical_front.Q_FYS1,
    &eqns_.p_.vertical_front.Q_FYS1);
  p_def_d(
    "tire.mf62_vertical.front.Q_FYS2", eqns_.p_.vertical_front.Q_FYS2,
    &eqns_.p_.vertical_front.Q_FYS2);
  p_def_d(
    "tire.mf62_vertical.front.Q_FYS3", eqns_.p_.vertical_front.Q_FYS3,
    &eqns_.p_.vertical_front.Q_FYS3);
  p_def_d(
    "tire.mf62_vertical.front.PFZ1", eqns_.p_.vertical_front.PFZ1, &eqns_.p_.vertical_front.PFZ1);
  p_def_d(
    "tire.mf62_vertical.front.BREFF", eqns_.p_.vertical_front.BREFF,
    &eqns_.p_.vertical_front.BREFF);
  p_def_d(
    "tire.mf62_vertical.front.DREFF", eqns_.p_.vertical_front.DREFF,
    &eqns_.p_.vertical_front.DREFF);
  p_def_d(
    "tire.mf62_vertical.front.FREFF", eqns_.p_.vertical_front.FREFF,
    &eqns_.p_.vertical_front.FREFF);
  p_def_d(
    "tire.mf62_vertical.front.BOTTOM_OFFST", eqns_.p_.vertical_front.BOTTOM_OFFST,
    &eqns_.p_.vertical_front.BOTTOM_OFFST);
  p_def_d(
    "tire.mf62_vertical.front.BOTTOM_STIFF", eqns_.p_.vertical_front.BOTTOM_STIFF,
    &eqns_.p_.vertical_front.BOTTOM_STIFF);
  p_def_d(
    "tire.mf62_vertical.front.QSY1", eqns_.p_.vertical_front.QSY1, &eqns_.p_.vertical_front.QSY1);
  p_def_d(
    "tire.mf62_vertical.front.QSY2", eqns_.p_.vertical_front.QSY2, &eqns_.p_.vertical_front.QSY2);
  p_def_d(
    "tire.mf62_vertical.front.QSY3", eqns_.p_.vertical_front.QSY3, &eqns_.p_.vertical_front.QSY3);
  p_def_d(
    "tire.mf62_vertical.front.QSY4", eqns_.p_.vertical_front.QSY4, &eqns_.p_.vertical_front.QSY4);
  p_def_d(
    "tire.mf62_vertical.front.QSY5", eqns_.p_.vertical_front.QSY5, &eqns_.p_.vertical_front.QSY5);
  p_def_d(
    "tire.mf62_vertical.front.QSY6", eqns_.p_.vertical_front.QSY6, &eqns_.p_.vertical_front.QSY6);
  p_def_d(
    "tire.mf62_vertical.front.QSY7", eqns_.p_.vertical_front.QSY7, &eqns_.p_.vertical_front.QSY7);
  p_def_d(
    "tire.mf62_vertical.front.QSY8", eqns_.p_.vertical_front.QSY8, &eqns_.p_.vertical_front.QSY8);
  p_def_d(
    "tire.mf62_vertical.rear.UNLOADED_RADIUS", eqns_.p_.vertical_rear.UNLOADED_RADIUS,
    &eqns_.p_.vertical_rear.UNLOADED_RADIUS);
  p_def_d(
    "tire.mf62_vertical.rear.WIDTH", eqns_.p_.vertical_rear.WIDTH, &eqns_.p_.vertical_rear.WIDTH);
  p_def_d(
    "tire.mf62_vertical.rear.ASPECT_RATIO", eqns_.p_.vertical_rear.ASPECT_RATIO,
    &eqns_.p_.vertical_rear.ASPECT_RATIO);
  p_def_d(
    "tire.mf62_vertical.rear.RIM_RADIUS", eqns_.p_.vertical_rear.RIM_RADIUS,
    &eqns_.p_.vertical_rear.RIM_RADIUS);
  p_def_d(
    "tire.mf62_vertical.rear.LONGVL", eqns_.p_.vertical_rear.LONGVL,
    &eqns_.p_.vertical_rear.LONGVL);
  p_def_d(
    "tire.mf62_vertical.rear.FNOMIN", eqns_.p_.vertical_rear.FNOMIN,
    &eqns_.p_.vertical_rear.FNOMIN);
  p_def_d(
    "tire.mf62_vertical.rear.NOMPRES", eqns_.p_.vertical_rear.NOMPRES,
    &eqns_.p_.vertical_rear.NOMPRES);
  p_def_d(
    "tire.mf62_vertical.rear.VERTICAL_STIFFNESS", eqns_.p_.vertical_rear.VERTICAL_STIFFNESS,
    &eqns_.p_.vertical_rear.VERTICAL_STIFFNESS);
  p_def_d(
    "tire.mf62_vertical.rear.VERTICAL_DAMPING", eqns_.p_.vertical_rear.VERTICAL_DAMPING,
    &eqns_.p_.vertical_rear.VERTICAL_DAMPING);
  p_def_d(
    "tire.mf62_vertical.rear.Q_RE0", eqns_.p_.vertical_rear.Q_RE0, &eqns_.p_.vertical_rear.Q_RE0);
  p_def_d(
    "tire.mf62_vertical.rear.Q_V1", eqns_.p_.vertical_rear.Q_V1, &eqns_.p_.vertical_rear.Q_V1);
  p_def_d(
    "tire.mf62_vertical.rear.Q_V2", eqns_.p_.vertical_rear.Q_V2, &eqns_.p_.vertical_rear.Q_V2);
  p_def_d(
    "tire.mf62_vertical.rear.Q_FZ2", eqns_.p_.vertical_rear.Q_FZ2, &eqns_.p_.vertical_rear.Q_FZ2);
  p_def_d(
    "tire.mf62_vertical.rear.Q_FCX", eqns_.p_.vertical_rear.Q_FCX, &eqns_.p_.vertical_rear.Q_FCX);
  p_def_d(
    "tire.mf62_vertical.rear.Q_FCY", eqns_.p_.vertical_rear.Q_FCY, &eqns_.p_.vertical_rear.Q_FCY);
  p_def_d(
    "tire.mf62_vertical.rear.Q_FCY2", eqns_.p_.vertical_rear.Q_FCY2,
    &eqns_.p_.vertical_rear.Q_FCY2);
  p_def_d(
    "tire.mf62_vertical.rear.Q_CAM1", eqns_.p_.vertical_rear.Q_CAM1,
    &eqns_.p_.vertical_rear.Q_CAM1);
  p_def_d(
    "tire.mf62_vertical.rear.Q_CAM2", eqns_.p_.vertical_rear.Q_CAM2,
    &eqns_.p_.vertical_rear.Q_CAM2);
  p_def_d(
    "tire.mf62_vertical.rear.Q_CAM3", eqns_.p_.vertical_rear.Q_CAM3,
    &eqns_.p_.vertical_rear.Q_CAM3);
  p_def_d(
    "tire.mf62_vertical.rear.Q_FYS1", eqns_.p_.vertical_rear.Q_FYS1,
    &eqns_.p_.vertical_rear.Q_FYS1);
  p_def_d(
    "tire.mf62_vertical.rear.Q_FYS2", eqns_.p_.vertical_rear.Q_FYS2,
    &eqns_.p_.vertical_rear.Q_FYS2);
  p_def_d(
    "tire.mf62_vertical.rear.Q_FYS3", eqns_.p_.vertical_rear.Q_FYS3,
    &eqns_.p_.vertical_rear.Q_FYS3);
  p_def_d(
    "tire.mf62_vertical.rear.PFZ1", eqns_.p_.vertical_rear.PFZ1, &eqns_.p_.vertical_rear.PFZ1);
  p_def_d(
    "tire.mf62_vertical.rear.BREFF", eqns_.p_.vertical_rear.BREFF, &eqns_.p_.vertical_rear.BREFF);
  p_def_d(
    "tire.mf62_vertical.rear.DREFF", eqns_.p_.vertical_rear.DREFF, &eqns_.p_.vertical_rear.DREFF);
  p_def_d(
    "tire.mf62_vertical.rear.FREFF", eqns_.p_.vertical_rear.FREFF, &eqns_.p_.vertical_rear.FREFF);
  p_def_d(
    "tire.mf62_vertical.rear.BOTTOM_OFFST", eqns_.p_.vertical_rear.BOTTOM_OFFST,
    &eqns_.p_.vertical_rear.BOTTOM_OFFST);
  p_def_d(
    "tire.mf62_vertical.rear.BOTTOM_STIFF", eqns_.p_.vertical_rear.BOTTOM_STIFF,
    &eqns_.p_.vertical_rear.BOTTOM_STIFF);
  p_def_d(
    "tire.mf62_vertical.rear.QSY1", eqns_.p_.vertical_rear.QSY1, &eqns_.p_.vertical_rear.QSY1);
  p_def_d(
    "tire.mf62_vertical.rear.QSY2", eqns_.p_.vertical_rear.QSY2, &eqns_.p_.vertical_rear.QSY2);
  p_def_d(
    "tire.mf62_vertical.rear.QSY3", eqns_.p_.vertical_rear.QSY3, &eqns_.p_.vertical_rear.QSY3);
  p_def_d(
    "tire.mf62_vertical.rear.QSY4", eqns_.p_.vertical_rear.QSY4, &eqns_.p_.vertical_rear.QSY4);
  p_def_d(
    "tire.mf62_vertical.rear.QSY5", eqns_.p_.vertical_rear.QSY5, &eqns_.p_.vertical_rear.QSY5);
  p_def_d(
    "tire.mf62_vertical.rear.QSY6", eqns_.p_.vertical_rear.QSY6, &eqns_.p_.vertical_rear.QSY6);
  p_def_d(
    "tire.mf62_vertical.rear.QSY7", eqns_.p_.vertical_rear.QSY7, &eqns_.p_.vertical_rear.QSY7);
  p_def_d(
    "tire.mf62_vertical.rear.QSY8", eqns_.p_.vertical_rear.QSY8, &eqns_.p_.vertical_rear.QSY8);
  p_def_d(
    "tire.inflation_reference_absolute_Pa.front_left", 1.0,
    &eqns_.p_.inflation_reference_absolute_Pa.front_left);
  p_def_d(
    "tire.inflation_reference_absolute_Pa.front_right", 1.0,
    &eqns_.p_.inflation_reference_absolute_Pa.front_right);
  p_def_d(
    "tire.inflation_reference_absolute_Pa.rear_left", 1.0,
    &eqns_.p_.inflation_reference_absolute_Pa.rear_left);
  p_def_d(
    "tire.inflation_reference_absolute_Pa.rear_right", 1.0,
    &eqns_.p_.inflation_reference_absolute_Pa.rear_right);
  p_def_d("suspension.kinematics.enabled", 1.0, &eqns_.p_.kinematics_enabled);
  p_def_d(
    "suspension.kinematics.material_velocity_enabled", 1.0, &eqns_.p_.material_velocity_enabled,
    "Reciprocal contact material velocity for coupled suspension; zero is a historical diagnostic");
  p_def_d("suspension.corner_weights.enabled", 1.0, &eqns_.p_.corner_weights_enabled);
  p_def_d(
    "suspension.corner_weights.kg.front_left", eqns_.p_.corner_weights_kg.front_left,
    &eqns_.p_.corner_weights_kg.front_left);
  p_def_d(
    "suspension.corner_weights.kg.front_right", eqns_.p_.corner_weights_kg.front_right,
    &eqns_.p_.corner_weights_kg.front_right);
  p_def_d(
    "suspension.corner_weights.kg.rear_left", eqns_.p_.corner_weights_kg.rear_left,
    &eqns_.p_.corner_weights_kg.rear_left);
  p_def_d(
    "suspension.corner_weights.kg.rear_right", eqns_.p_.corner_weights_kg.rear_right,
    &eqns_.p_.corner_weights_kg.rear_right);
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.travel_m", &eqns_.p_.kinematic_travel,
    eqns_.p_.kinematic_travel, tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Coupled rigid-link geometry knots");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.steering_mean_angle_rad",
    &eqns_.p_.kinematic_mean_angle, eqns_.p_.kinematic_mean_angle,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Coupled rigid-link geometry knots");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.rack_m", &eqns_.p_.kinematic_rack,
    eqns_.p_.kinematic_rack, tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Coupled rigid-link geometry knots");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.control_m.front",
    &eqns_.p_.kinematic_front.control, eqns_.p_.kinematic_front.control,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.spring.front",
    &eqns_.p_.kinematic_front.spring, eqns_.p_.kinematic_front.spring,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.damper.front",
    &eqns_.p_.kinematic_front.damper, eqns_.p_.kinematic_front.damper,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.normal_x.front",
    &eqns_.p_.kinematic_front.normal_x, eqns_.p_.kinematic_front.normal_x,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.normal_y.front",
    &eqns_.p_.kinematic_front.normal_y, eqns_.p_.kinematic_front.normal_y,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.normal_z.front",
    &eqns_.p_.kinematic_front.normal_z, eqns_.p_.kinematic_front.normal_z,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.contact_x.front",
    &eqns_.p_.kinematic_front.contact_x, eqns_.p_.kinematic_front.contact_x,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.contact_y.front",
    &eqns_.p_.kinematic_front.contact_y, eqns_.p_.kinematic_front.contact_y,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.contact_z.front",
    &eqns_.p_.kinematic_front.contact_z, eqns_.p_.kinematic_front.contact_z,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.jac_x.front",
    &eqns_.p_.kinematic_front.jac_x, eqns_.p_.kinematic_front.jac_x,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.jac_y.front",
    &eqns_.p_.kinematic_front.jac_y, eqns_.p_.kinematic_front.jac_y,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.jac_z.front",
    &eqns_.p_.kinematic_front.jac_z, eqns_.p_.kinematic_front.jac_z,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.angjac_x.front",
    &eqns_.p_.kinematic_front.angjac_x, eqns_.p_.kinematic_front.angjac_x,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.angjac_y.front",
    &eqns_.p_.kinematic_front.angjac_y, eqns_.p_.kinematic_front.angjac_y,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.angjac_z.front",
    &eqns_.p_.kinematic_front.angjac_z, eqns_.p_.kinematic_front.angjac_z,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.rocker.front",
    &eqns_.p_.kinematic_front.rocker, eqns_.p_.kinematic_front.rocker,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.control_m.rear",
    &eqns_.p_.kinematic_rear.control, eqns_.p_.kinematic_rear.control,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.spring.rear",
    &eqns_.p_.kinematic_rear.spring, eqns_.p_.kinematic_rear.spring,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.damper.rear",
    &eqns_.p_.kinematic_rear.damper, eqns_.p_.kinematic_rear.damper,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.normal_x.rear",
    &eqns_.p_.kinematic_rear.normal_x, eqns_.p_.kinematic_rear.normal_x,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.normal_y.rear",
    &eqns_.p_.kinematic_rear.normal_y, eqns_.p_.kinematic_rear.normal_y,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.normal_z.rear",
    &eqns_.p_.kinematic_rear.normal_z, eqns_.p_.kinematic_rear.normal_z,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.contact_x.rear",
    &eqns_.p_.kinematic_rear.contact_x, eqns_.p_.kinematic_rear.contact_x,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.contact_y.rear",
    &eqns_.p_.kinematic_rear.contact_y, eqns_.p_.kinematic_rear.contact_y,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.contact_z.rear",
    &eqns_.p_.kinematic_rear.contact_z, eqns_.p_.kinematic_rear.contact_z,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.jac_x.rear",
    &eqns_.p_.kinematic_rear.jac_x, eqns_.p_.kinematic_rear.jac_x,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.jac_y.rear",
    &eqns_.p_.kinematic_rear.jac_y, eqns_.p_.kinematic_rear.jac_y,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.jac_z.rear",
    &eqns_.p_.kinematic_rear.jac_z, eqns_.p_.kinematic_rear.jac_z,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.angjac_x.rear",
    &eqns_.p_.kinematic_rear.angjac_x, eqns_.p_.kinematic_rear.angjac_x,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.angjac_y.rear",
    &eqns_.p_.kinematic_rear.angjac_y, eqns_.p_.kinematic_rear.angjac_y,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.angjac_z.rear",
    &eqns_.p_.kinematic_rear.angjac_z, eqns_.p_.kinematic_rear.angjac_z,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.rocker.rear",
    &eqns_.p_.kinematic_rear.rocker, eqns_.p_.kinematic_rear.rocker,
    tam::pmg::ParameterType::DOUBLE_ARRAY, "Selected DummyVehicle rigid-link map, rack-major order");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.arb.front", &eqns_.p_.kinematic_front.arb,
    eqns_.p_.kinematic_front.arb, tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Closed rocker/drop-link/bar mechanism normalized to documented ground stiffness");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.kinematics.arb.rear", &eqns_.p_.kinematic_rear.arb,
    eqns_.p_.kinematic_rear.arb, tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Closed rocker/drop-link/bar mechanism normalized to documented ground stiffness");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.heave.rear_bar_travel_m",
    &eqns_.p_.rear_third_bar_travel, eqns_.p_.rear_third_bar_travel,
    tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Closed rear third-slider motion from mean bar rotation");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.suspension.heave.rear_device_compression_m",
    &eqns_.p_.rear_third_device_compression, eqns_.p_.rear_third_device_compression,
    tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Closed rear third-slider motion from mean bar rotation");
  // Anti-Roll bar
  p_def_d(
    "suspension.antirollbar_virtual_spring_stiffness_Npm.front", 12000, &eqns_.p_.c_ar_f,
    "Considerung the travel at the virtual suspension point at 0.5 track width");  // [N/m]
  p_def_d(
    "suspension.antirollbar_virtual_spring_stiffness_Npm.rear", 5000, &eqns_.p_.c_ar_r,
    "Considerung the travel at the virtual suspension point at 0.5 track width");  // [N/m]

  // Add toe
  p_def_d("suspension.toe_out_rad.front_left", 0.0, &eqns_.p_.toe_out_rad.front_left);
  p_def_d("suspension.toe_out_rad.front_right", 0.0, &eqns_.p_.toe_out_rad.front_right);
  p_def_d("suspension.toe_out_rad.rear_left", 0.0, &eqns_.p_.toe_out_rad.rear_left);
  p_def_d("suspension.toe_out_rad.rear_right", 0.0, &eqns_.p_.toe_out_rad.rear_right);

  p_def_d("suspension.camber_out_rad.front_left", 0.0, &eqns_.p_.camber_out_rad.front_left);
  p_def_d("suspension.camber_out_rad.front_right", 0.0, &eqns_.p_.camber_out_rad.front_right);
  p_def_d("suspension.camber_out_rad.rear_left", 0.0, &eqns_.p_.camber_out_rad.rear_left);
  p_def_d("suspension.camber_out_rad.rear_right", 0.0, &eqns_.p_.camber_out_rad.rear_right);

  p_def_d("tire.contact.belt.enabled", 1.0, &eqns_.p_.belt_enabled);
  p_def_d("tire.contact.belt.mass_kg.front", 1.0, &eqns_.p_.belt_mass_f);
  p_def_d("tire.contact.belt.mass_kg.rear", 1.0, &eqns_.p_.belt_mass_r);
  p_def_d("tire.contact.belt.frequency_hz.front", 1.0, &eqns_.p_.belt_frequency_f);
  p_def_d("tire.contact.belt.frequency_hz.rear", 1.0, &eqns_.p_.belt_frequency_r);
  p_def_d("suspension.resonant_damper.enabled", 1.0, &eqns_.p_.resonant_enabled);
  p_def_d("tire.contact.belt.resonant.enabled", 1.0, &eqns_.p_.resonant_belt_enabled);
  p_def_d("tire.contact.belt.resonant.frequency_hz", 1.0, &eqns_.p_.resonant_belt_hz);
  p_def_d("suspension.resonant_damper.frequency_hz", 1.0, &eqns_.p_.resonant_hz);
  p_def_d("suspension.resonant_damper.bandwidth_hz", 1.0, &eqns_.p_.resonant_width_hz);
  p_def_d("suspension.resonant_damper.residual_fraction", 1.0, &eqns_.p_.resonant_residual);
  p_def_d("tire.contact.belt.maxwell.enabled", 1.0, &eqns_.p_.maxwell_enabled);
  p_def_d("tire.contact.belt.maxwell.frequency_hz", 1.0, &eqns_.p_.maxwell_hz);
  p_def_d("tire.contact.belt.maxwell.stiffness_fraction", 1.0, &eqns_.p_.maxwell_fraction);
  p_def_d("tire.contact.belt.maxwell.viscous_fraction", 1.0, &eqns_.p_.maxwell_viscous_fraction);
  p_def_d("tire.effective_radius_only.enabled", 1.0, &eqns_.p_.effective_radius_only);
  p_def_d("tire.longitudinal_relaxation_length_m", 1.0, &eqns_.p_.longitudinal_relaxation_length_m);
  p_def_d("tire.contact.belt.damping_ratio.front", 1.0, &eqns_.p_.belt_damping_ratio_f);
  p_def_d("tire.contact.belt.damping_ratio.rear", 1.0, &eqns_.p_.belt_damping_ratio_r);
  p_def_d("tire.contact.enabled", 1.0, &eqns_.p_.contact_enabled);
  p_def_d(
    "tire.contact.strict_surface_coverage", 1.0, &eqns_.p_.strict_surface_coverage,
    "Reject missing collider coverage instead of a discontinuous ribbon fallback");
  p_def_d("tire.contact.damping_Nspm.front", 1.0, &eqns_.p_.tire_damping_f);
  p_def_d("tire.contact.damping_Nspm.rear", 1.0, &eqns_.p_.tire_damping_r);
  p_def_d("tire.contact.width_m.front", 1.0, &eqns_.p_.tire_width_f);
  p_def_d("tire.contact.width_m.rear", 1.0, &eqns_.p_.tire_width_r);
  p_def_d("tire.contact.half_length_m", 1.0, &eqns_.p_.contact_half_length);
  p_def_d("tire.contact.enveloping_radius_m.front", 1.0, &eqns_.p_.enveloping_radius_f);
  p_def_d("tire.contact.enveloping_radius_m.rear", 1.0, &eqns_.p_.enveloping_radius_r);
  p_def_d("tire.contact.micro_radius_m", 1.0, &eqns_.p_.contact_micro_radius);
  p_def_d("tire.contact.lateral_stiffness_Npm.front", 1.0, &eqns_.p_.lateral_stiffness_f);
  p_def_d("tire.contact.lateral_stiffness_Npm.rear", 1.0, &eqns_.p_.lateral_stiffness_r);
  p_def_d("tire.contact.sidewall_height_m.front", 1.0, &eqns_.p_.sidewall_height_f);
  p_def_d("tire.contact.sidewall_height_m.rear", 1.0, &eqns_.p_.sidewall_height_r);
  p_def_d(
    "tire.contact.lateral_drop_scale", 1.0, &eqns_.p_.lateral_drop_scale,
    "Geometric sidewall shortening; vertical/lateral coupling is a bounded hypothesis");
  p_def_d(
    "tire.contact.height_coupling_deflection_limit_m", 1.0, &eqns_.p_.lateral_drop_deflection_limit,
    "Bound the unvalidated vertical coupling at the tire-sheet 20-mm lateral range");
  p_def_d("suspension.bump_stop_gap_m.front", 1.0, &eqns_.p_.bump_gap_f);
  p_def_d("suspension.bump_stop_gap_m.rear", 1.0, &eqns_.p_.bump_gap_r);
  p_def_d("suspension.bump_stop_stiffness_Npm2.front", 1.0, &eqns_.p_.bump_stiffness_f);
  p_def_d("suspension.bump_stop_stiffness_Npm2.rear", 1.0, &eqns_.p_.bump_stiffness_r);

  p_def_d("suspension.heave.spring_stiffness_Npm.front", 1.0, &eqns_.p_.heave_spring_f);
  p_def_d("suspension.heave.spring_stiffness_Npm.rear", 1.0, &eqns_.p_.heave_spring_r);
  p_def_d("suspension.heave.rubber_stiffness_Npm2.front", 1.0, &eqns_.p_.heave_stiffness_f);
  p_def_d("suspension.heave.rubber_stiffness_Npm2.rear", 1.0, &eqns_.p_.heave_stiffness_r);
  p_def_d("suspension.heave.device_gap_m.front", 1.0, &eqns_.p_.heave_gap_f);
  p_def_d("suspension.heave.device_gap_m.rear", 1.0, &eqns_.p_.heave_gap_r);
  p_def_d("suspension.heave.wheel_device_motion_ratio.front", 1.0, &eqns_.p_.heave_motion_ratio_f);
  p_def_d("suspension.heave.wheel_device_motion_ratio.rear", 1.0, &eqns_.p_.heave_motion_ratio_r);

  // Tire characteristics
  p_def_d("tire.rolling_radius_m.front", 0.3, &eqns_.p_.rr_w_f, "At 20 mps");  // [m]@20mps
  p_def_d("tire.rolling_radius_m.rear", 0.3, &eqns_.p_.rr_w_r, "At 20 mps");   // [m]@20mps
  p_def_d("tire.spring_stiffness_Npm.front", 250000.0, &eqns_.p_.c_t_f);       // [N/m]
  p_def_d("tire.spring_stiffness_Npm.rear", 250000.0, &eqns_.p_.c_t_r);        // [N/m]
  p_def_d(
    "tire.rolling_resistance_coefficient", 0.025, &eqns_.p_.c_rr,
    "Rolling resistance coefficient. Ratio of F_roll to F_n");

  // Tire Radii
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.tire.rolling_radius_m.velocity_scaling.velocity_points_mps",
    &eqns_.p_.rr_vel_scale_vel_points_mps, std::vector<double>({0, 100}),
    tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Velocity points for scaling the tire radius dynamically with velocity");
  param_manager_->declare_parameter(
    "vehicle_dynamics_double_track.tire.rolling_radius_m.velocity_scaling.factors",
    &eqns_.p_.rr_vel_scale_scale_factors, std::vector<double>({1.0, 1.0}),
    tam::pmg::ParameterType::DOUBLE_ARRAY,
    "Scale factors for scaling the tire radius dynamically with velocity");

  // Parameters for Tire Model
  eqns_.tire_models_.front_left.declare_parameters(
    param_manager_.get(), "vehicle_dynamics_double_track.z.tire_model.front_left.");
  eqns_.tire_models_.front_right.declare_parameters(
    param_manager_.get(), "vehicle_dynamics_double_track.z.tire_model.front_right.");
  eqns_.tire_models_.rear_left.declare_parameters(
    param_manager_.get(), "vehicle_dynamics_double_track.z.tire_model.rear_left.");
  eqns_.tire_models_.rear_right.declare_parameters(
    param_manager_.get(), "vehicle_dynamics_double_track.z.tire_model.rear_right.");

  eqns_.aero_model_.declare_parameters(
    param_manager_.get(), "vehicle_dynamics_double_track.aerodynamics.");
}
}  // namespace tam::ocd::vehicle_dynamics
