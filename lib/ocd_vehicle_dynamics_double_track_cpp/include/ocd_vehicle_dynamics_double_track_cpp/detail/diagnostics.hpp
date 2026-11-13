// Copyright 2026 TUM Autonomous Motorsport
#pragma once
namespace tam::ocd::vehicle_dynamics
{
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_suspension_output()
{
  vd_output_.body_motion = {
    x_vec_[x::z_m],       x_vec_[x::v_z_mps],         x_dot_vec_[x::v_z_mps],
    x_vec_[x::phi_rad],   x_vec_[x::phi_dot_radps],   x_dot_vec_[x::phi_dot_radps],
    x_vec_[x::theta_rad], x_vec_[x::theta_dot_radps], x_dot_vec_[x::theta_dot_radps]};
  auto & out = vd_output_.suspension;
  out.wheel_displacement_m = {
    x_vec_[x::z_w_FL_m], x_vec_[x::z_w_FR_m], x_vec_[x::z_w_RL_m], x_vec_[x::z_w_RR_m]};
  // Recover the body attachment positions from the existing spring geometry.
  out.body_displacement_m = imr_.suspension_spring_compression_m +
                            imr_.suspension_spring_initial_compression_m +
                            out.wheel_displacement_m + kinematic_jacking_m_;
  if (p_.kinematics_enabled > .5)
    for (int i = 0; i < 4; ++i) {
      const double side = i % 2 ? -1. : 1., xx = i < 2 ? p_dep_.l_f_sprung : -p_dep_.l_r_sprung,
                   yy = side * (i < 2 ? p_.b_f : p_.b_r) / 2. - body_cg_y_m_;
      out.body_displacement_m[i] =
        x_vec_[x::z_m] + x_vec_[x::phi_rad] * yy - x_vec_[x::theta_rad] * xx;
    }
  out.body_displacement_over_road_m = out.body_displacement_m - imr_.road_contact_height_m;
  out.wheel_displacement_over_road_m = out.wheel_displacement_m - imr_.road_contact_height_m;
  out.wheel_velocity_mps = {
    x_vec_[x::v_z_w_FL_mps], x_vec_[x::v_z_w_FR_mps], x_vec_[x::v_z_w_RL_mps],
    x_vec_[x::v_z_w_RR_mps]};
  out.wheel_acceleration_mps2 = {
    x_dot_vec_[x::v_z_w_FL_mps], x_dot_vec_[x::v_z_w_FR_mps], x_dot_vec_[x::v_z_w_RL_mps],
    x_dot_vec_[x::v_z_w_RR_mps]};
  for (int i = 0; i < 4; ++i)
    out.damper_extension_m[i] = geometry_equivalent_extension_m_[i] /
                                (i < 2 ? p_.damper_motion_ratio_f : p_.damper_motion_ratio_r);
  if (p_.kinematics_enabled > .5) {
    const double f = -.5 * (kinematic_spring_m_[0] + kinematic_spring_m_[1]),
                 r = -p_.heave_motion_ratio_r * suspension_smooth_lookup(
                                                  .5 * (kinematic_arb_m_[2] + kinematic_arb_m_[3]),
                                                  p_.rear_third_bar_travel,
                                                  p_.rear_third_device_compression)
                                                  .first;
    out.third_equivalent_extension_m = {f, f, r, r};
  }
  out.spring_compression_m = imr_.suspension_spring_compression_m;
  out.damper_compression_velocity_mps = imr_.suspension_damper_compression_speed_mps;
  out.spring_force_N = imr_.suspension_spring_force_N;
  out.damper_force_N = imr_.suspension_damper_force_N;
  out.bump_force_N = imr_.suspension_bump_force_N + imr_.suspension_heave_force_N;
}
}  // namespace tam::ocd::vehicle_dynamics
