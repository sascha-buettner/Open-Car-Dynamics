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
  const auto positive = [](double value) { return std::isfinite(value) && value > 0.; };
  if (!(positive(p_.longitudinal_relaxation_length_m) &&
        (p_.longitudinal_recovery_length_front_m == 0. ||
         positive(p_.longitudinal_recovery_length_front_m)) &&
        (p_.longitudinal_recovery_length_rear_m == 0. ||
         positive(p_.longitudinal_recovery_length_rear_m)) &&
        ((p_.longitudinal_recovery_length_front_m == 0. &&
          p_.longitudinal_recovery_length_rear_m == 0.) ||
         positive(p_.longitudinal_recovery_gate_width)) &&
        positive(p_.longitudinal_relaxation_speed_min_mps) &&
        positive(p_.longitudinal_relaxation_speed_max_mps) &&
        p_.longitudinal_relaxation_speed_max_mps >= p_.longitudinal_relaxation_speed_min_mps &&
        positive(p_.lateral_relaxation_length_m) &&
        positive(p_.lateral_relaxation_speed_min_mps) &&
        positive(p_.lateral_relaxation_speed_max_mps) &&
        p_.lateral_relaxation_speed_max_mps >= p_.lateral_relaxation_speed_min_mps &&
        (p_.fy_relaxation_strength == 0. ||
         (positive(p_.lateral_relaxation_reference_front_N) &&
          positive(p_.lateral_relaxation_reference_rear_N))) &&
        (p_.sidewall_enabled < .5 ||
         (positive(p_.lateral_deflection_tau_min_s) &&
          positive(p_.lateral_deflection_length_m) &&
          positive(p_.lateral_deflection_speed_min_mps)))))
    throw std::invalid_argument("Explicit positive finite tire relaxation parameters required");
  auto pt1_T_relax_long = [this](double vx, double target, double state, bool front) {
    const double recovery_length = front ? p_.longitudinal_recovery_length_front_m :
                                           p_.longitudinal_recovery_length_rear_m;
    double length = p_.longitudinal_relaxation_length_m;
    if (recovery_length > 0. && state < 0. && target > state) {
      double gate = std::clamp(-state / p_.longitudinal_recovery_gate_width, 0., 1.);
      gate = gate * gate * (3. - 2. * gate);
      length += gate * (recovery_length - length);
    }
    return length /
           std::clamp(std::abs(vx), p_.longitudinal_relaxation_speed_min_mps,
                      p_.longitudinal_relaxation_speed_max_mps);
  };
  auto pt1_T_relax_lat = [this](double vx) {
    return p_.lateral_relaxation_length_m /
           std::clamp(std::abs(vx), p_.lateral_relaxation_speed_min_mps,
                      p_.lateral_relaxation_speed_max_mps);
  };

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
  if(p_.ring_enabled>.5){
    const auto a=ring_planar_acceleration();
    x_dot_vec_[x::v_x_mps]=a[0]+x_vec_[x::psi_dot_radps]*x_vec_[x::v_y_mps];
    x_dot_vec_[x::v_y_mps]=a[1]-x_vec_[x::psi_dot_radps]*x_vec_[x::v_x_mps];
    x_dot_vec_[x::psi_dot_radps]=a[2];
    if(p_.I_x<=4*p_.ring_inertia_trans||p_.I_y<=4*p_.ring_inertia_spin)throw std::domain_error("mfSCHWIFT chassis angular inertia partition invalid");
    x_dot_vec_[x::phi_dot_radps]=imr_.resulting_torque_Nm.x/(p_.I_x-4*p_.ring_inertia_trans);
    x_dot_vec_[x::theta_dot_radps]=imr_.resulting_torque_Nm.y/(p_.I_y-4*p_.ring_inertia_spin);
  }
  // Tire delay
  x_dot_vec_[x::kappa_tire_FL] = (imr_.tire_longitudinal_slip.front_left  - x_vec_[x::kappa_tire_FL]) / pt1_T_relax_long(imr_.velocity_wheel_over_ground_tire_frame_mps.front_left[0], imr_.tire_longitudinal_slip.front_left, x_vec_[x::kappa_tire_FL], true); // NOLINT
  x_dot_vec_[x::kappa_tire_FR] = (imr_.tire_longitudinal_slip.front_right - x_vec_[x::kappa_tire_FR]) / pt1_T_relax_long(imr_.velocity_wheel_over_ground_tire_frame_mps.front_right[0], imr_.tire_longitudinal_slip.front_right, x_vec_[x::kappa_tire_FR], true); // NOLINT
  x_dot_vec_[x::kappa_tire_RL] = (imr_.tire_longitudinal_slip.rear_left   - x_vec_[x::kappa_tire_RL]) / pt1_T_relax_long(imr_.velocity_wheel_over_ground_tire_frame_mps.rear_left[0], imr_.tire_longitudinal_slip.rear_left, x_vec_[x::kappa_tire_RL], false); // NOLINT
  x_dot_vec_[x::kappa_tire_RR] = (imr_.tire_longitudinal_slip.rear_right  - x_vec_[x::kappa_tire_RR]) / pt1_T_relax_long(imr_.velocity_wheel_over_ground_tire_frame_mps.rear_right[0], imr_.tire_longitudinal_slip.rear_right, x_vec_[x::kappa_tire_RR], false); // NOLINT
  x_dot_vec_[x::alpha_tire_FL] = (imr_.tire_slip_angle_rad.front_left     - x_vec_[x::alpha_tire_FL]) / (pt1_T_relax_lat(imr_.velocity_wheel_over_ground_tire_frame_mps.front_left[0])*lateral_relaxation_factor(imr_.contact_normal_load_N.front_left,p_.lateral_relaxation_reference_front_N,p_.fy_relaxation_strength)); // NOLINT
  x_dot_vec_[x::alpha_tire_FR] = (imr_.tire_slip_angle_rad.front_right    - x_vec_[x::alpha_tire_FR]) / (pt1_T_relax_lat(imr_.velocity_wheel_over_ground_tire_frame_mps.front_right[0])*lateral_relaxation_factor(imr_.contact_normal_load_N.front_right,p_.lateral_relaxation_reference_front_N,p_.fy_relaxation_strength)); // NOLINT
  x_dot_vec_[x::alpha_tire_RL] = (imr_.tire_slip_angle_rad.rear_left      - x_vec_[x::alpha_tire_RL]) / (pt1_T_relax_lat(imr_.velocity_wheel_over_ground_tire_frame_mps.rear_left[0])*lateral_relaxation_factor(imr_.contact_normal_load_N.rear_left,p_.lateral_relaxation_reference_rear_N,p_.fy_relaxation_strength)); // NOLINT
  x_dot_vec_[x::alpha_tire_RR] = (imr_.tire_slip_angle_rad.rear_right     - x_vec_[x::alpha_tire_RR]) / (pt1_T_relax_lat(imr_.velocity_wheel_over_ground_tire_frame_mps.rear_right[0])*lateral_relaxation_factor(imr_.contact_normal_load_N.rear_right,p_.lateral_relaxation_reference_rear_N,p_.fy_relaxation_strength)); // NOLINT
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
    if (p_.sidewall_enabled < .5) {
      x_dot_vec_[latidx[i]] = 0.;
      continue;
    }
    const double tau = std::max(
      p_.lateral_deflection_tau_min_s,
      p_.lateral_deflection_length_m /
        std::max(p_.lateral_deflection_speed_min_mps,
                 std::abs(imr_.velocity_wheel_over_ground_tire_frame_mps[i][0])));
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
                                       x_vec_[x::contact_capacity_FL_N + i],
                                       p_.load_memory_tau_s, p_.load_memory_release_tau_s)
                                   : 0.;
  if(p_.ring_enabled>.5 && p_.ring_residual_relaxation>.5)for(int i=0;i<4;++i){
    // Preserve supplied total quasi-static compliance: residual compliance is
    // the total minus all ring compliance seen at the ground contact point.
    const double mb=i<2?p_.belt_mass_f:p_.belt_mass_r,R=imr_.dynamic_tire_radius_m[i];
    const double kx=mb*std::pow(2*M_PI*(i<2?p_.belt_frequency_f:p_.belt_frequency_r),2);
    const double ky=mb*std::pow(2*M_PI*p_.ring_frequency_lat,2);
    const double kt=p_.ring_inertia_spin*std::pow(2*M_PI*p_.ring_frequency_windup,2);
    const double kg=p_.ring_inertia_trans*std::pow(2*M_PI*p_.ring_frequency_yaw,2);
    const double kcx=1./(1./kx+R*R/kt);
    const double kcy=1./(1./ky+std::pow((1.-p_.ring_lateral_shear_fraction)*R,2)/kg);
    const double ktx=p_.ring_longitudinal_stiffness,kty=i<2?p_.lateral_stiffness_f:p_.lateral_stiffness_r;
    if(!(kcx>ktx&&kcy>kty))throw std::domain_error("mfSCHWIFT nonpositive tangential residual compliance");
    const double krx=1./(1./ktx-1./kcx),kry=1./(1./kty-1./kcy);
    const double load=std::max(p_.low_load_boundary_N,imr_.contact_normal_load_N[i]);
    const double gamma=-(p_.camber_out_rad[i]+kinematic_camber_rad_[i]),pres=p_.inflation_absolute_Pa[i]-p_.atmosphere_Pa,eps=1e-5;
    const auto a=tire_models_[i].evaluate({eps,0.,load,gamma,pres});
    const auto b=tire_models_[i].evaluate({-eps,0.,load,gamma,pres});
    const auto c=tire_models_[i].evaluate({0.,eps,load,gamma,pres});
    const auto d=tire_models_[i].evaluate({0.,-eps,load,gamma,pres});
    const double Cx=std::abs((a.longitudinal_force_N-b.longitudinal_force_N)/(2*eps));
    const double Cy=std::abs((c.lateral_force_N-d.lateral_force_N)/(2*eps));
    const double speed=std::max(3.,std::abs(imr_.velocity_wheel_over_ground_tire_frame_mps[i][0]));
    // Below the measured load boundary the continuation has no identified
    // relaxation law. A declared 0.5-ms lower time bound regularizes this limit.
    const double tx=std::max(.0005,Cx/(krx*speed)),ty=std::max(.0005,Cy/(kry*speed));
    x_dot_vec_[x::kappa_tire_FL+i]=(imr_.tire_longitudinal_slip[i]-x_vec_[x::kappa_tire_FL+i])/tx;
    x_dot_vec_[x::alpha_tire_FL+i]=(imr_.tire_slip_angle_rad[i]-x_vec_[x::alpha_tire_FL+i])/ty;
  }
  const int ring_states[10]={x::ring_x_FL_m,x::ring_y_FL_m,x::ring_vx_FL_mps,x::ring_vy_FL_mps,x::ring_camber_FL_rad,x::ring_windup_FL_rad,x::ring_yaw_FL_rad,x::ring_omega_x_FL_radps,x::ring_omega_y_FL_radps,x::ring_omega_z_FL_radps};
  for(int j=0;j<10;++j)for(int i=0;i<4;++i)x_dot_vec_[ring_states[j]+i]=p_.ring_enabled>.5?ring_derivatives_[j][i]:0.;
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
