// Copyright 2026 TUM Autonomous Motorsport
#pragma once
namespace tam::ocd::vehicle_dynamics
{
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_dependent_parameters()
{
  p_dep_.l_r = p_.l - p_.l_f;
  p_dep_.m_sprung = p_.m - 2 * p_.m_t_f - 2 * p_.m_t_r;
  // Calculate the cog for the unsprung mass
  p_dep_.l_r_sprung = (p_.m * p_dep_.l_r - 2 * p_.m_t_f * p_.l) / p_dep_.m_sprung;
  p_dep_.l_f_sprung = p_.l - p_dep_.l_r_sprung;
  body_cg_y_m_ = vehicle_cg_y_m_ = 0.;
  if (p_.corner_weights_enabled > .5) {
    double total = 0., moment = 0., sprung_moment = 0.;
    for (int i = 0; i < 4; ++i) {
      const double side = i % 2 ? -1. : 1., y = side * (i < 2 ? p_.b_f : p_.b_r) / 2.;
      total += p_.corner_weights_kg[i];
      moment += p_.corner_weights_kg[i] * y;
      sprung_moment += (p_.corner_weights_kg[i] - (i < 2 ? p_.m_t_f : p_.m_t_r)) * y;
    }
    if (std::abs(total - p_.m) > 1e-5)
      throw std::invalid_argument("Corner weights must sum to physical mass");
    vehicle_cg_y_m_ = moment / p_.m;
    body_cg_y_m_ = sprung_moment / p_dep_.m_sprung;
  }
  // Calculate the parameters used for squat force calc
  if (p_.geometry_enabled > .5) {
    const double front_jounce = .5 * (x_vec_[x::z_w_FL_m] + x_vec_[x::z_w_FR_m]) - x_vec_[x::z_m] +
                                x_vec_[x::theta_rad] * p_dep_.l_f_sprung;
    const double rear_jounce = .5 * (x_vec_[x::z_w_RL_m] + x_vec_[x::z_w_RR_m]) - x_vec_[x::z_m] -
                               x_vec_[x::theta_rad] * p_dep_.l_r_sprung;
    p_.h_rc_front =
      suspension_lookup(front_jounce, p_.geometry_travel, p_.geometry_roll_front).first;
    p_.h_rc_rear = suspension_lookup(rear_jounce, p_.geometry_travel, p_.geometry_roll_rear).first;
  }
  p_dep_.slf =
    (((p_.geometry_enabled > .5 || p_.static_roll_centers_enabled > .5) ? p_.h_rc_front : p_.h_rc) /
     (0.5 * p_.b_f));
  p_dep_.slr =
    (((p_.geometry_enabled > .5 || p_.static_roll_centers_enabled > .5) ? p_.h_rc_rear : p_.h_rc) /
     (0.5 * p_.b_r));
  p_dep_.sadf_deceleration = (p_.h_pc_decel / (p_.l_f));
  p_dep_.sadr_deceleration = (p_.h_pc_decel / (p_dep_.l_r));
  p_dep_.sadf_acceleration = (p_.h_pc_accel / (p_.l_f));
  p_dep_.sadr_acceleration = (p_.h_pc_accel / (p_dep_.l_r));
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_effective_steering_angle()
{
  if (p_.kinematics_enabled > .5) {
    const double mean =
      .5 * (steering_angle_per_wheel_rad_.front_left + steering_angle_per_wheel_rad_.front_right);
    const auto rack = suspension_lookup(mean, p_.kinematic_mean_angle, p_.kinematic_rack);
    kinematic_rack_m_ = rack.first;
    kinematic_rack_rate_mps_ = rack.second * steering_mean_rate_radps_;
    const int zz[4] = {x::z_w_FL_m, x::z_w_FR_m, x::z_w_RL_m, x::z_w_RR_m};
    const int vv[4] = {x::v_z_w_FL_mps, x::v_z_w_FR_mps, x::v_z_w_RL_mps, x::v_z_w_RR_mps};
    const Eigen::Matrix3d body =
      (Eigen::AngleAxisd(x_vec_[x::theta_rad], Eigen::Vector3d::UnitY()) *
       Eigen::AngleAxisd(x_vec_[x::phi_rad], Eigen::Vector3d::UnitX()))
        .toRotationMatrix();
    for (int i = 0; i < 4; ++i) {
      const double side = i % 2 ? -1. : 1., xx = i < 2 ? p_dep_.l_f_sprung : -p_dep_.l_r_sprung,
                   yy = side * (i < 2 ? p_.b_f : p_.b_r) / 2. - body_cg_y_m_;
      const double attachment =
        x_vec_[x::z_m] + x_vec_[x::phi_rad] * yy - x_vec_[x::theta_rad] * xx;
      const double q = x_vec_[zz[i]] - attachment, control = i < 2 ? -side * kinematic_rack_m_ : 0.;
      const auto & table = i < 2 ? p_.kinematic_front : p_.kinematic_rear;
      auto lookup = [&](const auto & field) {
        return suspension_map(q, control, p_.kinematic_travel, table.control, field);
      };
      const auto spring = lookup(table.spring), damper = lookup(table.damper),
                 arb = lookup(table.arb);
      kinematic_arb_m_[i] = arb.value;
      kinematic_arb_jac_[i] = arb.dq;
      kinematic_q_m_[i] = q;
      kinematic_q_velocity_mps_[i] = x_vec_[vv[i]] - x_vec_[x::v_z_mps] -
                                     x_vec_[x::phi_dot_radps] * yy +
                                     x_vec_[x::theta_dot_radps] * xx;
      kinematic_spring_m_[i] = spring.value;
      kinematic_spring_jac_[i] = spring.dq;
      kinematic_spring_steer_jac_[i] = spring.ds;
      geometry_equivalent_extension_m_[i] = -damper.value;
      geometry_leverage_[i] = damper.dq;
      kinematic_damper_steer_jac_[i] = damper.ds;
      if (!(spring.dq > 0. && damper.dq > 0.))
        throw std::invalid_argument("Selected suspension mechanism reverses leverage");
      Eigen::Vector3d n(
        lookup(table.normal_x).value, lookup(table.normal_y).value, lookup(table.normal_z).value);
      Eigen::Vector3d cp(
        lookup(table.contact_x).value, lookup(table.contact_y).value,
        lookup(table.contact_z).value);
      Eigen::Vector3d j(
        lookup(table.jac_x).value, lookup(table.jac_y).value, lookup(table.jac_z).value);
      Eigen::Vector3d w(
        lookup(table.angjac_x).value, lookup(table.angjac_y).value, lookup(table.angjac_z).value);
      if (i % 2 == 0) {
        n.x() *= -1.;
        n.z() *= -1.;
        cp.y() *= -1.;
        j.y() *= -1.;
        w.x() *= -1.;
        w.z() *= -1.;
      }
      n = body *
          (Eigen::AngleAxisd(side * p_.toe_out_rad[i], Eigen::Vector3d::UnitZ()) * n.normalized());
      n.normalize();
      imr_.effective_steering_angle_per_wheel_rad[i] = std::atan2(-n.x(), n.y());
      const double gamma = -side * std::asin(std::clamp(n.z(), -1., 1.));
      kinematic_camber_rad_[i] = gamma - p_.camber_out_rad[i];
      tire_gamma_rad_[i] = -gamma;
      kinematic_toe_rad_[i] = imr_.effective_steering_angle_per_wheel_rad[i] - (i < 2 ? mean : 0.) -
                              side * p_.toe_out_rad[i];
      j = body * j;
      w = body * w;
      cp = body * cp;
      kinematic_contact_x_[i] = cp.x();
      kinematic_contact_y_[i] = cp.y();
      kinematic_contact_z_[i] = cp.z();
      kinematic_jac_x_[i] = j.x();
      kinematic_jac_y_[i] = j.y();
      kinematic_jac_z_[i] = j.z();
      kinematic_angjac_x_[i] = w.x();
      kinematic_angjac_y_[i] = w.y();
      kinematic_angjac_z_[i] = w.z();
    }
    for (int axle = 0; axle < 2; ++axle) {
      const int l = 2 * axle, r = l + 1;
      const double sl = kinematic_jac_y_[l] / kinematic_jac_z_[l],
                   sr = kinematic_jac_y_[r] / kinematic_jac_z_[r];
      const double yl = kinematic_contact_y_[l], yr = kinematic_contact_y_[r],
                   zl = external_influences_.z_height_road_m[l],
                   zr = external_influences_.z_height_road_m[r];
      const double cy = std::abs(sl - sr) > 1e-9 ? (zl - zr + sl * yl - sr * yr) / (sl - sr) : 0.;
      const double cz = zl - sl * (cy - yl) - .5 * (zl + zr);
      if (axle == 0) {
        roll_center_front_y_m_ = cy;
        p_.h_rc_front = cz;
      } else {
        roll_center_rear_y_m_ = cy;
        p_.h_rc_rear = cz;
      }
    }
    return;
  }
  double_per_wheel_t toe_addition = p_.toe_out_rad;
  if (p_.geometry_enabled > .5) {
    const double wheelz[4] = {
      x_vec_[x::z_w_FL_m], x_vec_[x::z_w_FR_m], x_vec_[x::z_w_RL_m], x_vec_[x::z_w_RR_m]};
    for (int i = 0; i < 4; ++i) {
      const double side = i % 2 == 0 ? 1. : -1.;
      const double attachment =
        x_vec_[x::z_m] + side * x_vec_[x::phi_rad] * (i < 2 ? p_.b_f : p_.b_r) / 2. +
        x_vec_[x::theta_rad] * (i < 2 ? -p_dep_.l_f_sprung : p_dep_.l_r_sprung);
      const double jounce = wheelz[i] - attachment;
      kinematic_toe_rad_[i] =
        suspension_lookup(
          jounce, p_.geometry_travel, i < 2 ? p_.geometry_toe_front : p_.geometry_toe_rear)
          .first;
      kinematic_camber_rad_[i] =
        suspension_lookup(
          jounce, p_.geometry_travel, i < 2 ? p_.geometry_camber_front : p_.geometry_camber_rear)
          .first;
      // Orient the wheel normal about the tilted steering axis (Rodrigues).
      // Camber is measured relative to the ribbon, including sprung-body roll.
      kinematic_camber_rad_[i] -= side * x_vec_[x::phi_rad];
      if (i < 2) {
        Eigen::Vector3d axis(
          -std::tan(p_.caster_rad), -side * std::tan(p_.kingpin_inclination_rad), 1.);
        axis.normalize();
        const double delta = steering_angle_per_wheel_rad_[i];
        const double nz = axis.x() * std::sin(delta) + axis.z() * axis.y() * (1. - std::cos(delta));
        kinematic_camber_rad_[i] -= side * std::asin(std::clamp(nz, -1., 1.));
      }
      toe_addition[i] += kinematic_toe_rad_[i];
    }
  }
  // Toe out means steering to right on the right side -> Therefore negate the sign.
  toe_addition.front_right *= -1;
  toe_addition.rear_right *= -1;

  imr_.effective_steering_angle_per_wheel_rad = steering_angle_per_wheel_rad_ + toe_addition;
}
}  // namespace tam::ocd::vehicle_dynamics
