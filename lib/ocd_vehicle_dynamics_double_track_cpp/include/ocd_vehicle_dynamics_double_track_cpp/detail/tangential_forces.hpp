// Copyright 2026 TUM Autonomous Motorsport
#pragma once
namespace tam::ocd::vehicle_dynamics
{
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<
  TireModelT, AeroModelT>::calculate_velocity_wheel_over_ground_mps()
{  // Tire Velocity Vector @ x - y axis
  imr_.velocity_wheel_over_ground_mps.front_left = {
    x_vec_[x::v_x_mps] - p_.b_f * x_vec_[x::psi_dot_radps] / 2,
    x_vec_[x::v_y_mps] + p_.l_f * x_vec_[x::psi_dot_radps]};
  imr_.velocity_wheel_over_ground_mps.front_right = {
    x_vec_[x::v_x_mps] + p_.b_f * x_vec_[x::psi_dot_radps] / 2,
    x_vec_[x::v_y_mps] + p_.l_f * x_vec_[x::psi_dot_radps]};
  imr_.velocity_wheel_over_ground_mps.rear_left = {
    x_vec_[x::v_x_mps] - p_.b_r * x_vec_[x::psi_dot_radps] / 2,
    x_vec_[x::v_y_mps] - p_dep_.l_r * x_vec_[x::psi_dot_radps]};
  imr_.velocity_wheel_over_ground_mps.rear_right = {
    x_vec_[x::v_x_mps] + p_.b_r * x_vec_[x::psi_dot_radps] / 2,
    x_vec_[x::v_y_mps] - p_dep_.l_r * x_vec_[x::psi_dot_radps]};
  if (p_.kinematics_enabled > .5 && p_.material_velocity_enabled > .5) {
    // The force calculation uses J^T F, so slip must use the same material J v.
    // Wheel spin remains the separate relative rotor coordinate R*omega.
    const double h = p_.h_cg + x_vec_[x::z_m];
    for (int i = 0; i < 4; ++i) {
      const double xplan = i < 2 ? p_.l_f : -p_dep_.l_r;
      imr_.velocity_wheel_over_ground_mps[i] = {
        x_vec_[x::v_x_mps] -
          (kinematic_contact_y_[i] - vehicle_cg_y_m_) * x_vec_[x::psi_dot_radps] -
          h * x_vec_[x::theta_dot_radps] + kinematic_jac_x_[i] * kinematic_q_velocity_mps_[i],
        x_vec_[x::v_y_mps] + (xplan + kinematic_contact_x_[i]) * x_vec_[x::psi_dot_radps] +
          h * x_vec_[x::phi_dot_radps] + kinematic_jac_y_[i] * kinematic_q_velocity_mps_[i]};
    }
  }
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<
  TireModelT, AeroModelT>::calculate_velocity_wheel_over_ground_tire_frame_mps()
{
  // Rotate the wheel velocity over ground into the vehicle coordinate frame (passive rotation,
  // i.e., reverse rotation)

  // Required for macro below
  double angle = 0;
  Eigen::Matrix2d rot_matrix;

#define EVAL_VELOCITY_WHEEL_OVER_GROUND_TIRE_FRAME(wheel)        \
                                                                 \
  angle = imr_.effective_steering_angle_per_wheel_rad.wheel;     \
  rot_matrix << cos(angle), sin(angle), -sin(angle), cos(angle); \
  imr_.velocity_wheel_over_ground_tire_frame_mps.wheel =         \
    rot_matrix * imr_.velocity_wheel_over_ground_mps.wheel;

  EVAL_MACRO_PER_WHEEL(EVAL_VELOCITY_WHEEL_OVER_GROUND_TIRE_FRAME);
#undef EVAL_VELOCITY_WHEEL_OVER_GROUND_TIRE_FRAME
  const int wheel_z_rate[4] = {
    p_.belt_enabled > .5 ? x::v_z_belt_FL_mps : x::v_z_w_FL_mps,
    p_.belt_enabled > .5 ? x::v_z_belt_FR_mps : x::v_z_w_FR_mps,
    p_.belt_enabled > .5 ? x::v_z_belt_RL_mps : x::v_z_w_RL_mps,
    p_.belt_enabled > .5 ? x::v_z_belt_RR_mps : x::v_z_w_RR_mps};
  for (int i = 0; i < 4; ++i) {
    const auto horizontal = imr_.velocity_wheel_over_ground_tire_frame_mps[i];
    imr_.velocity_wheel_over_ground_tire_frame_mps[i] = contact_frames_[i].velocity(
      {horizontal[0], horizontal[1],
       x_vec_[wheel_z_rate[i]] + ((p_.kinematics_enabled > .5 && p_.material_velocity_enabled > .5)
                                    ? (kinematic_jac_z_[i] - 1.) * kinematic_q_velocity_mps_[i]
                                    : 0.)});
  }
  // clang-format on
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_velocity_tire_rotation_mps()
{
  // Velocity due to tire roll @ x - y axis:
  imr_.velocity_tire_rotation_mps.front_left =  // p_.rr_w_f * wheel_speeds_radps_.front_left;
    imr_.dynamic_tire_radius_m.front_left * wheel_speeds_radps_.front_left;
  imr_.velocity_tire_rotation_mps.front_right =  // p_.rr_w_f * wheel_speeds_radps_.front_right;
    imr_.dynamic_tire_radius_m.front_right * wheel_speeds_radps_.front_right;
  imr_.velocity_tire_rotation_mps.rear_left =  // p_.rr_w_r * wheel_speeds_radps_.rear_left;
    imr_.dynamic_tire_radius_m.rear_left * wheel_speeds_radps_.rear_left;
  imr_.velocity_tire_rotation_mps.rear_right =  // p_.rr_w_r * wheel_speeds_radps_.rear_right;
    imr_.dynamic_tire_radius_m.rear_right * wheel_speeds_radps_.rear_right;
  for (int i = 0; i < 4; ++i)
    imr_.velocity_tire_rotation_mps[i] *= contact_frames_[i].spin_lever_scale;
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_tire_longitudinal_slip()
{
  // TODO(Simon S) verify slip definition

  constexpr static auto long_slip_def =
    [](double rotational_velocity_mps, double vel_over_ground_mps) {
      return std::clamp(
        (rotational_velocity_mps - vel_over_ground_mps) / std::max(3.0, vel_over_ground_mps), -1.0,
        1.0);
    };

#define CALC_LONG_SLIP(x)                        \
  imr_.tire_longitudinal_slip.x = long_slip_def( \
    imr_.velocity_tire_rotation_mps.x, imr_.velocity_wheel_over_ground_tire_frame_mps.x[0]);

  EVAL_MACRO_PER_WHEEL(CALC_LONG_SLIP);
#undef CALC_LONG_SLIP
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_tire_slip_angle_rad()
{
  auto slip_angle_def = [this](double vy, double vx) { return -atan2(vy, std::max(1.0, vx)); };

  // clang-format off
  #define CALC_TIRE_SLIP(x)                               \
  imr_.tire_slip_angle_rad.x = slip_angle_def(           \
    imr_.velocity_wheel_over_ground_tire_frame_mps.x[1], \
    imr_.velocity_wheel_over_ground_tire_frame_mps.x[0])

  EVAL_MACRO_PER_WHEEL(CALC_TIRE_SLIP);
  #undef CALC_TIRE_SLIP
  // clang-format on
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_tire_forces_tire_frame_N()
{
  tire_mx_Nm_ = {0., 0., 0., 0.};
  tire_mz_Nm_ = {0., 0., 0., 0.};
  // Zero out tire forces when at rest
  if (
    p_.force_conditioning_enabled < .5 && std::pow(x_vec_[x::v_x_mps], 2) < 0.001 &&
    std::pow(x_vec_[x::v_y_mps], 2) < 0.001 &&
    std::pow(wheel_speeds_radps_.front_left, 2) < 0.001 &&
    std::pow(wheel_speeds_radps_.front_right, 2) < 0.001 &&
    std::pow(wheel_speeds_radps_.rear_left, 2) < 0.001 &&
    std::pow(wheel_speeds_radps_.rear_right, 2) < 0.001) {
    imr_.tire_forces_tire_frame_N.front_right = {0, 0};
    imr_.tire_forces_tire_frame_N.front_left = {0, 0};
    imr_.tire_forces_tire_frame_N.rear_left = {0, 0};
    imr_.tire_forces_tire_frame_N.rear_right = {0, 0};
    return;
  }

  auto tire_model_eval = [this](
                           double F_z_N, double long_slip, double slip_angle_rad,
                           TireModelT & model, Eigen::Vector2d & target, bool mirror,
                           double external_scale, double camber_out_rad, double pressure_Pa,
                           double & mx, double & mz, int wheel) {
    // Airborne tires transmit no friction force. MF52 divides by Fz and is
    // undefined at zero load; do not introduce a fictitious contact preload.
    if (F_z_N <= 0.0) {
      target.setZero();
      return;
    }
    // Mirror the tire model to account for asymmetry
    int mirror_factor = (mirror) ? -1 : 1;
    if (p_.force_conditioning_enabled > .5) {
      long_slip = p_.force_condition_kappa[wheel];
      slip_angle_rad = p_.force_condition_alpha[wheel];
    }
    const double force_load = std::max(F_z_N, p_.low_load_boundary_N);
    auto out = model.evaluate(
      {long_slip, mirror_factor * slip_angle_rad, force_load, -camber_out_rad, pressure_Pa});
    const double low_scale = F_z_N / force_load;
    const double memory =
      p_.load_memory_fraction > 0.
        ? load_capacity_factor(
            F_z_N, x_vec_[x::contact_capacity_FL_N + wheel], p_.load_memory_fraction)
        : 1.;
    out.longitudinal_force_N *= low_scale * memory;
    out.lateral_force_N *= low_scale * memory;
    out.self_aligning_moment_Nm *= low_scale * memory;
    out.overturning_moment_Nm *= low_scale;
    target[0] = out.longitudinal_force_N * external_scale;
    target[1] = mirror_factor * out.lateral_force_N * external_scale;
    mx = mirror_factor * out.overturning_moment_Nm * external_scale;
    mz = mirror_factor * out.self_aligning_moment_Nm * external_scale;
  };

  // Use the delayed kappa and alpha
  double_per_wheel_t kappa, alpha;
  kappa.front_left = x_vec_[x::kappa_tire_FL];
  kappa.front_right = x_vec_[x::kappa_tire_FR];
  kappa.rear_left = x_vec_[x::kappa_tire_RL];
  kappa.rear_right = x_vec_[x::kappa_tire_RR];
  alpha.front_left = x_vec_[x::alpha_tire_FL];
  alpha.front_right = x_vec_[x::alpha_tire_FR];
  alpha.rear_left = x_vec_[x::alpha_tire_RL];
  alpha.rear_right = x_vec_[x::alpha_tire_RR];

  // clang-format off
  #define EVAL_TIRE_MODEL(x, mirror, wheel_index)   \
  tire_model_eval( \
    imr_.contact_normal_load_N.x, \
    kappa.x, \
    alpha.x, \
    tire_models_.x, \
    imr_.tire_forces_tire_frame_N.x, \
    mirror, \
    external_influences_.lambda_mue.x, \
    (p_.camber_out_rad.x + kinematic_camber_rad_.x), (p_.inflation_absolute_Pa.x-p_.atmosphere_Pa),tire_mx_Nm_.x,tire_mz_Nm_.x,wheel_index \
  )

  EVAL_TIRE_MODEL(front_left, false, 0);
  EVAL_TIRE_MODEL(rear_left, false, 2);
  EVAL_TIRE_MODEL(front_right, true, 1);
  EVAL_TIRE_MODEL(rear_right, true, 3);
  #undef EVAL_TIRE_MODEL
  // clang-format on
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_tire_forces_N()
{
  for (int i = 0; i < 4; ++i) {
    const auto tangent = contact_frames_[i].force(
      imr_.tire_forces_tire_frame_N[i][0], imr_.tire_forces_tire_frame_N[i][1]);
    const double delta = imr_.effective_steering_angle_per_wheel_rad[i];
    imr_.tire_forces_N[i] = {
      std::cos(delta) * tangent.x() - std::sin(delta) * tangent.y(),
      std::sin(delta) * tangent.x() + std::cos(delta) * tangent.y()};
    imr_.contact_tangential_vertical_force_N[i] = tangent.z();
  }
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_tire_rolling_resistance_N()
{
  if (p_.mf62_vertical_enabled > .5 && p_.mf62_rolling_enabled > .5) {
    for (int i = 0; i < 4; ++i) {
      const auto & t = i < 2 ? p_.vertical_front : p_.vertical_rear;
      imr_.tire_rolling_resistance_N[i] =
        mf62_rolling_torque(
          t, imr_.vertical_tire_force_N[i], p_.inflation_absolute_Pa[i] - p_.atmosphere_Pa,
          x_vec_[x::v_x_mps], tire_gamma_rad_[i], imr_.tire_forces_tire_frame_N[i][0]) /
        imr_.dynamic_tire_radius_m[i];
    }
    return;
  }

// Always positive since always going forward
#define CALC_TIRE_RR(x_in) \
  imr_.tire_rolling_resistance_N.x_in = imr_.vertical_tire_force_N.x_in * p_.c_rr;
  EVAL_MACRO_PER_WHEEL(CALC_TIRE_RR);
#undef CALC_TIRE_RR
}
}  // namespace tam::ocd::vehicle_dynamics
