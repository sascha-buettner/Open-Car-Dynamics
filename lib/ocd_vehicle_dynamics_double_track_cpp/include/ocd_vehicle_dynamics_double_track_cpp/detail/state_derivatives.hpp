// Copyright 2026 TUM Autonomous Motorsport
#pragma once
namespace tam::ocd::vehicle_dynamics
{
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::evaluate()
{
  // Calc the required forces and moments
  calculate_intermediate_results();
  calculate_drivetrain_load_torque_Nm();
  calculate_steering_load_torque_Nm();

  // According to Pacejka 2006 page 343, Eq. (7.18)
  if (!(std::isfinite(p_.longitudinal_relaxation_length_m) &&
        p_.longitudinal_relaxation_length_m > 0.))
    throw std::invalid_argument("Positive finite longitudinal relaxation length required");
  auto pt1_T_relax_long = [this](double vx) {
    return p_.longitudinal_relaxation_length_m / std::clamp(std::abs(vx), 1.0, 1.0);
  };
  auto pt1_T_relax_lat = [](double vx) { return 1.0 / std::clamp(std::abs(vx), 1.0, 1.0); };

  // State Equation
  // ============================================================================================
  // clang-format off
  // Vehicle Movement
  x_dot_vec_[x::x_m] = cos(x_vec_[x::psi_rad]) * x_vec_[x::v_x_mps] - sin(x_vec_[x::psi_rad]) * x_vec_[x::v_y_mps]; // NOLINT
  x_dot_vec_[x::y_m] = sin(x_vec_[x::psi_rad]) * x_vec_[x::v_x_mps] + cos(x_vec_[x::psi_rad]) * x_vec_[x::v_y_mps];  // NOLINT
  x_dot_vec_[x::psi_rad] = x_vec_[x::psi_dot_radps];  // psi_dot
  x_dot_vec_[x::v_x_mps] = imr_.resulting_force_N.x / p_.m + x_vec_[x::psi_dot_radps] * x_vec_[x::v_y_mps];  // NOLINT
  x_dot_vec_[x::v_y_mps] = imr_.resulting_force_N.y / p_.m - x_vec_[x::psi_dot_radps] * x_vec_[x::v_x_mps]; // NOLINT
  x_dot_vec_[x::psi_dot_radps] = imr_.resulting_torque_Nm.z / p_.I_z;
  // Vertical movement of vehicle body
  x_dot_vec_[x::z_m] = x_vec_[x::v_z_mps];
  x_dot_vec_[x::v_z_mps] = imr_.resulting_force_N.z / p_dep_.m_sprung;
  // Pitching an rolling
  x_dot_vec_[x::phi_rad] = x_vec_[x::phi_dot_radps];
  x_dot_vec_[x::phi_dot_radps] = imr_.resulting_torque_Nm.x / p_.I_x;
  x_dot_vec_[x::theta_rad] = x_vec_[x::theta_dot_radps];
  x_dot_vec_[x::theta_dot_radps] = imr_.resulting_torque_Nm.y / p_.I_y;
  // Wheel vertical dynamics
  x_dot_vec_[x::z_w_FL_m] = x_vec_[x::v_z_w_FL_mps];  // z_t_dot_LF
  x_dot_vec_[x::z_w_FR_m] = x_vec_[x::v_z_w_FR_mps];  // z_t_dot_RF
  x_dot_vec_[x::z_w_RL_m] = x_vec_[x::v_z_w_RL_mps];  // z_t_dot_LR
  x_dot_vec_[x::z_w_RR_m] = x_vec_[x::v_z_w_RR_mps];  // z_t_dot_RR
  x_dot_vec_[x::v_z_w_FL_mps] = imr_.resulting_vertical_force_on_wheel_N.front_left / (p_.m_t_f-(p_.belt_enabled>.5?p_.belt_mass_f:0.));
  x_dot_vec_[x::v_z_w_FR_mps] = imr_.resulting_vertical_force_on_wheel_N.front_right / (p_.m_t_f-(p_.belt_enabled>.5?p_.belt_mass_f:0.));
  x_dot_vec_[x::v_z_w_RL_mps] = imr_.resulting_vertical_force_on_wheel_N.rear_left / (p_.m_t_r-(p_.belt_enabled>.5?p_.belt_mass_r:0.));
  x_dot_vec_[x::v_z_w_RR_mps] = imr_.resulting_vertical_force_on_wheel_N.rear_right / (p_.m_t_r-(p_.belt_enabled>.5?p_.belt_mass_r:0.));
  // Tire delay
  x_dot_vec_[x::kappa_tire_FL] = (imr_.tire_longitudinal_slip.front_left  - x_vec_[x::kappa_tire_FL]) / pt1_T_relax_long(imr_.velocity_wheel_over_ground_tire_frame_mps.front_left[0]); // NOLINT
  x_dot_vec_[x::kappa_tire_FR] = (imr_.tire_longitudinal_slip.front_right - x_vec_[x::kappa_tire_FR]) / pt1_T_relax_long(imr_.velocity_wheel_over_ground_tire_frame_mps.front_right[0]); // NOLINT
  x_dot_vec_[x::kappa_tire_RL] = (imr_.tire_longitudinal_slip.rear_left   - x_vec_[x::kappa_tire_RL]) / pt1_T_relax_long(imr_.velocity_wheel_over_ground_tire_frame_mps.rear_left[0]); // NOLINT
  x_dot_vec_[x::kappa_tire_RR] = (imr_.tire_longitudinal_slip.rear_right  - x_vec_[x::kappa_tire_RR]) / pt1_T_relax_long(imr_.velocity_wheel_over_ground_tire_frame_mps.rear_right[0]); // NOLINT
  x_dot_vec_[x::alpha_tire_FL] = (imr_.tire_slip_angle_rad.front_left     - x_vec_[x::alpha_tire_FL]) / (pt1_T_relax_lat(imr_.velocity_wheel_over_ground_tire_frame_mps.front_left[0])*lateral_relaxation_factor(imr_.contact_normal_load_N.front_left,1.0,p_.fy_relaxation_strength)); // NOLINT
  x_dot_vec_[x::alpha_tire_FR] = (imr_.tire_slip_angle_rad.front_right    - x_vec_[x::alpha_tire_FR]) / (pt1_T_relax_lat(imr_.velocity_wheel_over_ground_tire_frame_mps.front_right[0])*lateral_relaxation_factor(imr_.contact_normal_load_N.front_right,1.0,p_.fy_relaxation_strength)); // NOLINT
  x_dot_vec_[x::alpha_tire_RL] = (imr_.tire_slip_angle_rad.rear_left      - x_vec_[x::alpha_tire_RL]) / (pt1_T_relax_lat(imr_.velocity_wheel_over_ground_tire_frame_mps.rear_left[0])*lateral_relaxation_factor(imr_.contact_normal_load_N.rear_left,1.0,p_.fy_relaxation_strength)); // NOLINT
  x_dot_vec_[x::alpha_tire_RR] = (imr_.tire_slip_angle_rad.rear_right     - x_vec_[x::alpha_tire_RR]) / (pt1_T_relax_lat(imr_.velocity_wheel_over_ground_tire_frame_mps.rear_right[0])*lateral_relaxation_factor(imr_.contact_normal_load_N.rear_right,1.0,p_.fy_relaxation_strength)); // NOLINT
  // clang-format on

  const int belt_z[4] = {x::z_belt_FL_m, x::z_belt_FR_m, x::z_belt_RL_m, x::z_belt_RR_m};
  const int belt_v[4] = {
    x::v_z_belt_FL_mps, x::v_z_belt_FR_mps, x::v_z_belt_RL_mps, x::v_z_belt_RR_mps};
  for (int i = 0; i < 4; ++i) {
    x_dot_vec_[belt_z[i]] = p_.belt_enabled > .5 ? x_vec_[belt_v[i]] : 0.;
    x_dot_vec_[belt_v[i]] = p_.belt_enabled > .5 ? imr_.belt_resulting_force_N[i] /
                                                     (i < 2 ? p_.belt_mass_f : p_.belt_mass_r)
                                                 : 0.;
  }
  const int latidx[4] = {
    x::lateral_deflection_tire_FL_m, x::lateral_deflection_tire_FR_m,
    x::lateral_deflection_tire_RL_m, x::lateral_deflection_tire_RR_m};
  for (int i = 0; i < 4; ++i) {
    const double tau = std::max(
      1.0, 1.0 / std::max(1., std::abs(imr_.velocity_wheel_over_ground_tire_frame_mps[i][0])));
    x_dot_vec_[latidx[i]] = (imr_.tire_lateral_target_m[i] - x_vec_[latidx[i]]) / tau;
  }
  for (int i = 0; i < 4; ++i) {
    x_dot_vec_[x::resonant_q_FL + i] = (p_.resonant_enabled > .5 || p_.resonant_belt_enabled > .5)
                                         ? x_vec_[x::resonant_v_FL + i]
                                         : 0.;
    x_dot_vec_[x::resonant_v_FL + i] = (p_.resonant_enabled > .5 || p_.resonant_belt_enabled > .5)
                                         ? p_.resonant_acceleration[i]
                                         : 0.;
    x_dot_vec_[x::belt_maxwell_FL + i] = p_.maxwell_enabled > .5 ? p_.maxwell_derivative[i] : 0.;
  }
  for (int i = 0; i < 4; ++i)
    x_dot_vec_[x::contact_capacity_FL_N + i] =
      p_.load_memory_fraction > 0. ? load_capacity_derivative(
                                       imr_.contact_normal_load_N[i],
                                       x_vec_[x::contact_capacity_FL_N + i], p_.load_memory_tau_s)
                                   : 0.;
  // Output
  // ============================================================================================
  calculate_suspension_output();
  vd_output_.position_m.x = x_vec_[x::x_m];
  vd_output_.position_m.y = x_vec_[x::y_m];
  vd_output_.orientation_rad.z = x_vec_[x::psi_rad];
  vd_output_.velocity_mps.x = x_vec_[x::v_x_mps];
  vd_output_.velocity_mps.y = x_vec_[x::v_y_mps];
  vd_output_.angular_velocity_radps.z = x_vec_[x::psi_dot_radps];
  vd_output_.angular_acceleration_radps2.z = x_dot_vec_[x::psi_dot_radps];
  vd_output_.tire_longitudinal_slip_per_wheel.front_left = x_vec_[x::kappa_tire_FL];
  vd_output_.tire_longitudinal_slip_per_wheel.front_right = x_vec_[x::kappa_tire_FR];
  vd_output_.tire_longitudinal_slip_per_wheel.rear_left = x_vec_[x::kappa_tire_RL];
  vd_output_.tire_longitudinal_slip_per_wheel.rear_right = x_vec_[x::kappa_tire_RR];
  vd_output_.tire_slip_angle_per_wheel_rad.front_left = x_vec_[x::alpha_tire_FL];
  vd_output_.tire_slip_angle_per_wheel_rad.front_right = x_vec_[x::alpha_tire_FR];
  vd_output_.tire_slip_angle_per_wheel_rad.rear_left = x_vec_[x::alpha_tire_RL];
  vd_output_.tire_slip_angle_per_wheel_rad.rear_right = x_vec_[x::alpha_tire_RR];
#define GET_TIRE_FORCES(wheel)                                      \
  vd_output_.longitudinal_tire_force_tire_frame_per_wheel_N.wheel = \
    imr_.tire_forces_tire_frame_N.wheel[0];                         \
  vd_output_.lateral_tire_force_tire_frame_per_wheel_N.wheel =      \
    imr_.tire_forces_tire_frame_N.wheel[1];
  GET_TIRE_FORCES(front_left);
  GET_TIRE_FORCES(front_right);
  GET_TIRE_FORCES(rear_left);
  GET_TIRE_FORCES(rear_right);
#undef GET_TIRE_FORCES
  vd_output_.vertical_tire_force_per_wheel_N =
    imr_.vertical_tire_force_N + imr_.contact_tangential_vertical_force_N;
  //
  // Adapt the local accelerations to global accelerations since this is what the imu would measure
  vd_output_.acceleration_mps2.x =
    x_dot_vec_[x::v_x_mps] - x_vec_[x::psi_dot_radps] * x_vec_[x::v_y_mps];
  vd_output_.acceleration_mps2.y =
    x_dot_vec_[x::v_y_mps] + x_vec_[x::psi_dot_radps] * x_vec_[x::v_x_mps];
  vd_output_.acceleration_mps2.z = tam::constants::g_earth;
};
}  // namespace tam::ocd::vehicle_dynamics
