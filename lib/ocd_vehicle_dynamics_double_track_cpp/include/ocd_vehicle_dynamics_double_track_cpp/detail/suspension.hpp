// Copyright 2026 TUM Autonomous Motorsport
#pragma once
namespace tam::ocd::vehicle_dynamics
{
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<
  TireModelT, AeroModelT>::calculate_suspension_spring_initial_compression_m()
{
  if (p_.corner_weights_enabled > .5) {
    for (int i = 0; i < 4; ++i) {
      const auto & t = i < 2 ? p_.kinematic_front : p_.kinematic_rear;
      const double jz = p_.kinematics_enabled > .5
                          ? suspension_map(0., 0., p_.kinematic_travel, t.control, t.jac_z).value
                          : 1.;
      imr_.suspension_spring_initial_compression_m[i] =
        (p_.corner_weights_kg[i] * jz - (i < 2 ? p_.m_t_f : p_.m_t_r)) * G_CONST /
        (i < 2 ? p_.c_f : p_.c_r);
    }
    return;
  }
  // Compression of suspension coils @ standstill
  imr_.suspension_spring_initial_compression_m.front_left =
    p_dep_.m_sprung * G_CONST * p_dep_.l_r_sprung / (2 * p_.c_f * p_.l);
  imr_.suspension_spring_initial_compression_m.rear_left =
    p_dep_.m_sprung * G_CONST * p_dep_.l_f_sprung / (2 * p_.c_r * p_.l);
  imr_.suspension_spring_initial_compression_m.front_right =
    imr_.suspension_spring_initial_compression_m.front_left;
  imr_.suspension_spring_initial_compression_m.rear_right =
    imr_.suspension_spring_initial_compression_m.rear_left;
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<
  TireModelT, AeroModelT>::calculate_suspension_spring_compression_m()
{
  if (p_.kinematics_enabled > .5) {
    for (int i = 0; i < 4; ++i)
      imr_.suspension_spring_compression_m[i] =
        -imr_.suspension_spring_initial_compression_m[i] - kinematic_q_m_[i];
    return;
  }
  imr_.suspension_spring_compression_m.front_left =
    -imr_.suspension_spring_initial_compression_m.front_left + x_vec_[x::z_m] +
    x_vec_[x::phi_rad] * p_.b_f / 2 - x_vec_[x::theta_rad] * p_dep_.l_f_sprung -
    x_vec_[x::z_w_FL_m];
  imr_.suspension_spring_compression_m.front_right =
    -imr_.suspension_spring_initial_compression_m.front_right + x_vec_[x::z_m] -
    x_vec_[x::phi_rad] * p_.b_f / 2 - x_vec_[x::theta_rad] * p_dep_.l_f_sprung -
    x_vec_[x::z_w_FR_m];
  imr_.suspension_spring_compression_m.rear_left =
    -imr_.suspension_spring_initial_compression_m.rear_left + x_vec_[x::z_m] +
    x_vec_[x::phi_rad] * p_.b_r / 2 + x_vec_[x::theta_rad] * p_dep_.l_r_sprung -
    x_vec_[x::z_w_RL_m];
  imr_.suspension_spring_compression_m.rear_right =
    -imr_.suspension_spring_initial_compression_m.rear_right + x_vec_[x::z_m] -
    x_vec_[x::phi_rad] * p_.b_r / 2 + x_vec_[x::theta_rad] * p_dep_.l_r_sprung -
    x_vec_[x::z_w_RR_m];
  kinematic_jacking_m_ = {0., 0., 0., 0.};
  if (p_.geometry_enabled > .5 && p_.kinematics_enabled < .5) {
    const double mean =
      .5 * (steering_angle_per_wheel_rad_.front_left + steering_angle_per_wheel_rad_.front_right);
    for (int i = 0; i < 2; ++i) {
      const double signed_angle = (i == 0 ? -1. : 1.) * mean;
      const double jack =
        suspension_lookup(signed_angle, p_.geometry_steer, p_.geometry_jacking).first;
      kinematic_jacking_m_[i] = jack;
      // A body rise of +jack leaves the damper unchanged at fixed wheel height.
      imr_.suspension_spring_compression_m[i] -= jack;
    }
  }
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<
  TireModelT, AeroModelT>::calculate_suspension_damper_compression_speed_mps()
{
  if (p_.kinematics_enabled > .5) {
    for (int i = 0; i < 4; ++i) {
      const double rate = i < 2 ? (i % 2 ? 1. : -1.) * kinematic_rack_rate_mps_ : 0.;
      imr_.suspension_damper_compression_speed_mps[i] =
        -kinematic_q_velocity_mps_[i] -
        kinematic_damper_steer_jac_[i] * rate / geometry_leverage_[i];
    }
    return;
  }
  // Push/ Pull of suspension dampers
  imr_.suspension_damper_compression_speed_mps.front_left =
    x_vec_[x::v_z_mps] + x_vec_[x::phi_dot_radps] * p_.b_f / 2 -
    x_vec_[x::theta_dot_radps] * p_dep_.l_f_sprung - x_vec_[x::v_z_w_FL_mps];
  imr_.suspension_damper_compression_speed_mps.front_right =
    x_vec_[x::v_z_mps] - x_vec_[x::phi_dot_radps] * p_.b_f / 2 -
    x_vec_[x::theta_dot_radps] * p_dep_.l_f_sprung - x_vec_[x::v_z_w_FR_mps];
  imr_.suspension_damper_compression_speed_mps.rear_left =
    x_vec_[x::v_z_mps] + x_vec_[x::phi_dot_radps] * p_.b_r / 2 +
    x_vec_[x::theta_dot_radps] * p_dep_.l_r_sprung - x_vec_[x::v_z_w_RL_mps];
  imr_.suspension_damper_compression_speed_mps.rear_right =
    x_vec_[x::v_z_mps] - x_vec_[x::phi_dot_radps] * p_.b_r / 2 +
    x_vec_[x::theta_dot_radps] * p_dep_.l_r_sprung - x_vec_[x::v_z_w_RR_mps];
  if (p_.geometry_enabled > .5 && p_.kinematics_enabled < .5) {
    const double mean =
      .5 * (steering_angle_per_wheel_rad_.front_left + steering_angle_per_wheel_rad_.front_right);
    for (int i = 0; i < 2; ++i) {
      const double side = i == 0 ? -1. : 1.;
      imr_.suspension_damper_compression_speed_mps[i] -=
        side * steering_mean_rate_radps_ *
        suspension_lookup(side * mean, p_.geometry_steer, p_.geometry_jacking).second;
    }
  }
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_suspension_spring_force_N()
{
  if (p_.kinematics_enabled > .5) {
    for (int i = 0; i < 4; ++i) {
      const auto & t = i < 2 ? p_.kinematic_front : p_.kinematic_rear;
      const double j0 = suspension_map(0., 0., p_.kinematic_travel, t.control, t.spring).dq,
                   k = i < 2 ? p_.c_f : p_.c_r;
      const double preload = imr_.suspension_spring_initial_compression_m[i] / j0;
      imr_.suspension_spring_force_N[i] =
        -k * (kinematic_spring_m_[i] + preload) * kinematic_spring_jac_[i];
    }
    return;
  }

  // Suspension spring forces
  imr_.suspension_spring_force_N.front_left =
    p_.c_f * imr_.suspension_spring_compression_m.front_left;
  imr_.suspension_spring_force_N.front_right =
    p_.c_f * imr_.suspension_spring_compression_m.front_right;
  imr_.suspension_spring_force_N.rear_left =
    p_.c_r * imr_.suspension_spring_compression_m.rear_left;
  imr_.suspension_spring_force_N.rear_right =
    p_.c_r * imr_.suspension_spring_compression_m.rear_right;
  geometry_leverage_ = {1., 1., 1., 1.};
  geometry_equivalent_extension_m_ =
    imr_.suspension_spring_compression_m + imr_.suspension_spring_initial_compression_m;
  if (p_.geometry_enabled > .5)
    for (int i = 0; i < 4; ++i) {
      const double q = geometry_equivalent_extension_m_[i];
      auto g = suspension_smooth_lookup(
        -q, p_.geometry_travel, i < 2 ? p_.geometry_equiv_front : p_.geometry_equiv_rear);
      if (!(g.second > 0. && std::isfinite(g.first)))
        throw std::invalid_argument("Suspension lever must be positive");
      geometry_equivalent_extension_m_[i] = -g.first;
      geometry_leverage_[i] = g.second;
      const double k = i < 2 ? p_.c_f : p_.c_r;
      // Potential U=-F_static*q + .5*k*g(q)^2; preserve nominal equilibrium.
      imr_.suspension_spring_force_N[i] =
        -k * imr_.suspension_spring_initial_compression_m[i] + k * (-g.first) * g.second;
    }
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_suspension_damper_force_N()
{
  imr_.damper_curve_speed_extrapolated = {0, 0, 0, 0};
  for (int axle = 0; axle < 2; ++axle) {
    const bool enabled = (axle == 0 ? p_.damper_curve_enabled_f : p_.damper_curve_enabled_r) > .5;
    const double motion_ratio = axle == 0 ? p_.damper_motion_ratio_f : p_.damper_motion_ratio_r;
    const auto & speed = p_.damper_speed_mps;
    const auto & bump = axle == 0 ? p_.damper_bump_f : p_.damper_bump_r;
    const auto & rebound = axle == 0 ? p_.damper_rebound_f : p_.damper_rebound_r;
    if (enabled) {
      if (
        !(std::isfinite(motion_ratio) && motion_ratio > 0.) || speed.size() < 2 ||
        bump.size() != speed.size() || rebound.size() != speed.size() || speed.front() != 0. ||
        bump.front() != 0. || rebound.front() != 0.)
        throw std::invalid_argument("Damper curve dimensions/zero origin/motion ratio invalid");
      for (std::size_t i = 0; i < speed.size(); ++i)
        if (
          !(std::isfinite(speed[i]) && std::isfinite(bump[i]) && std::isfinite(rebound[i]) &&
            bump[i] >= 0. && rebound[i] >= 0.) ||
          (i && (speed[i] <= speed[i - 1] || bump[i] < bump[i - 1] || rebound[i] < rebound[i - 1])))
          throw std::invalid_argument("Damper curves must be finite, nonnegative and monotone");
    }
    for (int side = 0; side < 2; ++side) {
      const int i = 2 * axle + side;
      const double velocity =
        imr_.suspension_damper_compression_speed_mps[i] * geometry_leverage_[i];
      if (!enabled) {
        imr_.suspension_damper_force_N[i] =
          (axle == 0 ? p_.d_f : p_.d_r) * velocity * geometry_leverage_[i];
        continue;
      }
      // Positive native compression-speed coordinate means physical extension.
      // The same wheel/shaft ratio converts speed and reciprocal force; thus
      // F_wheel * v_wheel = F_shaft * v_shaft and damping is always dissipative.
      const double shaft_speed = std::abs(velocity) / motion_ratio;
      const auto & force = velocity >= 0. ? rebound : bump;
      const auto next = std::upper_bound(speed.begin(), speed.end(), shaft_speed);
      const auto hi = std::clamp<std::size_t>(next - speed.begin(), 1, speed.size() - 1);
      const double value = force[hi - 1] + (shaft_speed - speed[hi - 1]) *
                                             (force[hi] - force[hi - 1]) /
                                             (speed[hi] - speed[hi - 1]);
      imr_.suspension_damper_force_N[i] =
        std::copysign(value * geometry_leverage_[i] / motion_ratio, velocity);
      imr_.damper_curve_speed_extrapolated[i] = shaft_speed > speed.back() ? 1. : 0.;
      if (p_.resonant_enabled > .5) {
        if (!(p_.resonant_hz > 0. && p_.resonant_width_hz > 0. && p_.resonant_residual >= 0. &&
              p_.resonant_residual < 1.))
          throw std::invalid_argument("Invalid passive resonant damper");
        const int iq = x::resonant_q_FL + i, iv = x::resonant_v_FL + i;
        const double vr = x_vec_[iv], vrel = velocity - vr;
        const auto & curve = vrel >= 0. ? rebound : bump;
        const double ss = std::abs(vrel) / motion_ratio;
        const auto upper = std::upper_bound(speed.begin(), speed.end(), ss);
        const auto ih = std::clamp<std::size_t>(upper - speed.begin(), 1, speed.size() - 1);
        const double fv = curve[ih - 1] + (ss - speed[ih - 1]) * (curve[ih] - curve[ih - 1]) /
                                            (speed[ih] - speed[ih - 1]);
        const double branch = (1. - p_.resonant_residual) * std::copysign(fv / motion_ratio, vrel);
        const double c0 = .5 * (bump[1] + rebound[1]) / speed[1] / (motion_ratio * motion_ratio) *
                          (1. - p_.resonant_residual);
        const double inertia = c0 / (2 * M_PI * p_.resonant_width_hz),
                     stiffness = inertia * std::pow(2 * M_PI * p_.resonant_hz, 2);
        p_.resonant_acceleration[i] = (branch - stiffness * x_vec_[iq]) / inertia;
        imr_.suspension_damper_force_N[i] =
          p_.resonant_residual * imr_.suspension_damper_force_N[i] + branch * geometry_leverage_[i];
      }
    }
  }
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_antiroll_bar_force_N()
{
  if (p_.kinematics_enabled > .5) {
    for (int axle = 0; axle < 2; ++axle) {
      const int l = 2 * axle, r = l + 1;
      const double K = axle == 0 ? p_.c_ar_f : p_.c_ar_r,
                   twist = kinematic_arb_m_[l] - kinematic_arb_m_[r];
      imr_.antiroll_bar_force_N[l] = -.5 * K * twist * kinematic_arb_jac_[l];
      imr_.antiroll_bar_force_N[r] = .5 * K * twist * kinematic_arb_jac_[r];
    }
    return;
  }

  // Anti-Roll bar forces
  imr_.antiroll_bar_force_N.front_left = 0.5 * p_.c_ar_f *
                                         (imr_.suspension_spring_compression_m.front_left -
                                          imr_.suspension_spring_compression_m.front_right);
  imr_.antiroll_bar_force_N.rear_left = 0.5 * p_.c_ar_r *
                                        (imr_.suspension_spring_compression_m.rear_left -
                                         imr_.suspension_spring_compression_m.rear_right);
  imr_.antiroll_bar_force_N.front_right = -imr_.antiroll_bar_force_N.front_left;
  imr_.antiroll_bar_force_N.rear_right = -imr_.antiroll_bar_force_N.rear_left;
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_axle_vertical_force_N()
{  // Total Axle Forces
  if (p_.kinematics_enabled > .5) {
    for (int i = 0; i < 4; ++i)
      imr_.axle_vertical_force_N[i] =
        imr_.tire_forces_N[i][0] * kinematic_jac_x_[i] +
        imr_.tire_forces_N[i][1] * kinematic_jac_y_[i] +
        imr_.vertical_tire_force_N[i] * (kinematic_jac_z_[i] - 1.) +
        tire_mx_Nm_[i] *
          (std::cos(imr_.effective_steering_angle_per_wheel_rad[i]) * kinematic_angjac_x_[i] +
           std::sin(imr_.effective_steering_angle_per_wheel_rad[i]) * kinematic_angjac_y_[i]) +
        tire_mz_Nm_[i] * kinematic_angjac_z_[i];
    return;
  }

  // See page A-30
  // https://babel.hathitrust.org/cgi/pt?id=mdp.39015075298698&view=1up&seq=147&skin=2021
  auto calc_axle_force = [this](
                           double F_x, double F_y, double K_sli, double K_sadi, double K_sadi2) {
    double k_sadi_used = (F_x > 0.0) ? K_sadi2 : K_sadi;
    return (K_sli * F_y - k_sadi_used * F_x);
  };

  imr_.axle_vertical_force_N.front_left = calc_axle_force(
    -imr_.tire_forces_N.front_left[0], imr_.tire_forces_N.front_left[1], p_dep_.slf,
    p_dep_.sadf_deceleration, p_dep_.sadf_acceleration);
  imr_.axle_vertical_force_N.front_right = calc_axle_force(
    -imr_.tire_forces_N.front_right[0], -imr_.tire_forces_N.front_right[1], p_dep_.slf,
    p_dep_.sadf_deceleration, p_dep_.sadf_acceleration);
  imr_.axle_vertical_force_N.rear_left = calc_axle_force(
    imr_.tire_forces_N.rear_left[0], imr_.tire_forces_N.rear_left[1], p_dep_.slr,
    p_dep_.sadr_deceleration, p_dep_.sadr_acceleration);
  imr_.axle_vertical_force_N.rear_right = calc_axle_force(
    imr_.tire_forces_N.rear_right[0], -imr_.tire_forces_N.rear_right[1], p_dep_.slr,
    p_dep_.sadr_deceleration, p_dep_.sadr_acceleration);
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<
  TireModelT, AeroModelT>::calculate_resulting_suspension_force_N()
{
  // clang-format off
  #define VERTICAL_FORCES_SUSPENSION(tire_in)                                        \
    imr_.resulting_suspension_force_N.tire_in = \
    (  imr_.suspension_spring_force_N.tire_in \
    + imr_.suspension_damper_force_N.tire_in \
    + imr_.antiroll_bar_force_N.tire_in \
    + imr_.suspension_bump_force_N.tire_in \
    + imr_.suspension_heave_force_N.tire_in \
    + imr_.axle_vertical_force_N.tire_in \
    )

  EVAL_MACRO_PER_WHEEL(VERTICAL_FORCES_SUSPENSION);
  #undef VERTICAL_FORCES_SUSPENSION
  // clang-format on
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<
  TireModelT, AeroModelT>::calculate_resulting_vertical_force_on_wheel_N()
{
  const int rim_z[4] = {x::z_w_FL_m, x::z_w_FR_m, x::z_w_RL_m, x::z_w_RR_m};
  const int rim_v[4] = {x::v_z_w_FL_mps, x::v_z_w_FR_mps, x::v_z_w_RL_mps, x::v_z_w_RR_mps};
  const int belt_z[4] = {x::z_belt_FL_m, x::z_belt_FR_m, x::z_belt_RL_m, x::z_belt_RR_m};
  const int belt_v[4] = {
    x::v_z_belt_FL_mps, x::v_z_belt_FR_mps, x::v_z_belt_RL_mps, x::v_z_belt_RR_mps};
  for (int i = 0; i < 4; ++i) {
    const double total_mass = i < 2 ? p_.m_t_f : p_.m_t_r;
    double rim_mass = total_mass;
    double support = -imr_.tire_spring_force_N[i] + imr_.contact_tangential_vertical_force_N[i];
    imr_.belt_resulting_force_N[i] = 0.;
    imr_.belt_carcass_force_N[i] = 0.;
    imr_.belt_elastic_energy_J[i] = 0.;
    if (p_.belt_enabled > .5) {
      const double mb = i < 2 ? p_.belt_mass_f : p_.belt_mass_r,
                   frequency = i < 2 ? p_.belt_frequency_f : p_.belt_frequency_r;
      const double ratio = i < 2 ? p_.belt_damping_ratio_f : p_.belt_damping_ratio_r;
      const double cb = mb * std::pow(2. * M_PI * frequency, 2),
                   db = 2. * ratio * mb * 2. * M_PI * frequency;
      rim_mass -= mb;
      const double nominal_load =
        p_.radial_lookup_enabled > .5
          ? radial_preload_N_[i]
          : imr_.tire_spring_initial_compression_m[i] * (i < 2 ? p_.c_t_f : p_.c_t_r);
      const double compression = x_vec_[belt_z[i]] - x_vec_[rim_z[i]] +
                                 (nominal_load - mb * G_CONST) / cb - imr_.tire_lateral_drop_m[i];
      double carcass = cb * compression + db * (x_vec_[belt_v[i]] - x_vec_[rim_v[i]]);
      if (p_.maxwell_enabled > .5) {
        if (!(p_.maxwell_hz > 0. && p_.maxwell_fraction >= 0. && p_.maxwell_fraction <= .2 &&
              p_.maxwell_viscous_fraction >= 0. && p_.maxwell_viscous_fraction <= 1.))
          throw std::invalid_argument("Invalid Maxwell belt branch");
        const double dynamic_compression = x_vec_[belt_z[i]] - x_vec_[rim_z[i]];
        const double lag = x_vec_[x::belt_maxwell_FL + i];
        p_.maxwell_derivative[i] = 2 * M_PI * p_.maxwell_hz * (dynamic_compression - lag);
        carcass = cb * compression +
                  p_.maxwell_viscous_fraction * db * (x_vec_[belt_v[i]] - x_vec_[rim_v[i]]) +
                  p_.maxwell_fraction * cb * (dynamic_compression - lag);
      }
      if (p_.resonant_belt_enabled > .5) {
        if (p_.resonant_enabled > .5 || p_.maxwell_enabled > .5)
          throw std::invalid_argument("Use one memory branch per diagnostic");
        const double velocity = x_vec_[belt_v[i]] - x_vec_[rim_v[i]];
        const double inertia = (1. - p_.resonant_residual) * db / (2 * M_PI * p_.resonant_width_hz);
        const double stiffness = inertia * std::pow(2 * M_PI * p_.resonant_belt_hz, 2);
        const double branch =
          (1. - p_.resonant_residual) * db * (velocity - x_vec_[x::resonant_v_FL + i]);
        p_.resonant_acceleration[i] = (branch - stiffness * x_vec_[x::resonant_q_FL + i]) / inertia;
        carcass = cb * compression + p_.resonant_residual * db * velocity + branch;
      }
      imr_.belt_carcass_force_N[i] = carcass;
      imr_.belt_elastic_energy_J[i] = .5 * cb * compression * compression;
      imr_.belt_resulting_force_N[i] =
        support - carcass - mb * G_CONST + external_influences_.tire_force_z_N[i];
      support = carcass;
    }
    imr_.rim_support_force_N[i] = support;
    imr_.resulting_vertical_force_on_wheel_N[i] =
      imr_.resulting_suspension_force_N[i] + support - rim_mass * G_CONST +
      external_influences_.rim_force_z_N[i] +
      (p_.belt_enabled > .5 ? 0. : external_influences_.tire_force_z_N[i]);
  }
}
}  // namespace tam::ocd::vehicle_dynamics
