// Copyright 2026 TUM Autonomous Motorsport
#pragma once
namespace tam::ocd::vehicle_dynamics
{
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_aerodynamics()
{
  // Aerodynamic forces
  types::AeroModelInput input;
  input.vx_mps = x_vec_[x::v_x_mps];
  input.vy_mps = x_vec_[x::v_y_mps];
  input.wind_mps = external_influences_.wind_mps;
  input.pitch_angle_rad = x_vec_[x::theta_rad];
  input.z_m = x_vec_[x::z_m];
  input.pitch_rate_radps = x_vec_[x::theta_dot_radps];
  input.heave_velocity_mps = x_vec_[x::v_z_mps];

  types::AeroModelOutput aero_model_out = aero_model_.evaluate(input);
  imr_.aero_force_N = aero_model_out.force_cog_N;
  imr_.aero_torque_Nm = aero_model_out.torque_Nm;
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_resulting_force_N()
{
  // Scale external influences at low speed
  double external_influence_scale_factor_lowspeed =
    std::clamp((0.2 * (x_vec_[x::v_x_mps] - 2.0)), 0.0, 1.0);

  imr_.resulting_force_N.x =
    imr_.tire_forces_N.front_left[0] + imr_.tire_forces_N.front_right[0] +
    imr_.tire_forces_N.rear_left[0] + imr_.tire_forces_N.rear_right[0] + imr_.aero_force_N.x +
    external_influences_.external_force_N.x * external_influence_scale_factor_lowspeed;

  // Sum of tire forces @ y axis
  imr_.resulting_force_N.y =
    imr_.tire_forces_N.front_left[1] + imr_.tire_forces_N.front_right[1] +
    imr_.tire_forces_N.rear_left[1] + imr_.tire_forces_N.rear_right[1] + imr_.aero_force_N.y +
    external_influences_.external_force_N.y * external_influence_scale_factor_lowspeed;
  // Heave of Car
  // clang-format off
  imr_.resulting_force_N.z =
    - imr_.resulting_suspension_force_N.front_left
    - imr_.resulting_suspension_force_N.front_right
    - imr_.resulting_suspension_force_N.rear_left
    - imr_.resulting_suspension_force_N.rear_right
    - p_dep_.m_sprung * G_CONST
    + imr_.aero_force_N.z
    + external_influences_.external_force_N.z;
  // clang-format on
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_resulting_torque_Nm()
{
  // Sum of tire moments @ z axis
  // clang-format off
  imr_.resulting_torque_Nm.z =
      (imr_.tire_forces_N.front_left[1] + imr_.tire_forces_N.front_right[1]) * p_.l_f
    + (imr_.tire_forces_N.front_right[0] - imr_.tire_forces_N.front_left[0]) * p_.b_f / 2
    + (imr_.tire_forces_N.rear_right[0] - imr_.tire_forces_N.rear_left[0]) * p_.b_r / 2
    - (imr_.tire_forces_N.rear_left[1] + imr_.tire_forces_N.rear_right[1]) * p_dep_.l_r
    + imr_.aero_torque_Nm.z
    + external_influences_.external_torque_Nm.z;


  // Macro to sum of the forces of the spring, damper and antirollbar
  imr_.resulting_torque_Nm.y =
    - (imr_.tire_forces_N.front_left[0] + imr_.tire_forces_N.front_right[0] +
       imr_.tire_forces_N.rear_left[0] + imr_.tire_forces_N.rear_right[0]
      ) * p_.h_cg
    + imr_.resulting_suspension_force_N.front_left * p_.l_f
    + imr_.resulting_suspension_force_N.front_right * p_.l_f
    - imr_.resulting_suspension_force_N.rear_left * p_dep_.l_r
    - imr_.resulting_suspension_force_N.rear_right* p_dep_.l_r
    + p_dep_.m_sprung* G_CONST * (p_.l_f - p_dep_.l_f_sprung)
    + imr_.aero_torque_Nm.y
    + external_influences_.external_torque_Nm.y;

  // Sum of sprung mass moments @ x axis
  imr_.resulting_torque_Nm.x = (
       imr_.tire_forces_N.front_left[1] + imr_.tire_forces_N.front_right[1] +
       imr_.tire_forces_N.rear_left[1] + imr_.tire_forces_N.rear_right[1]
      ) * p_.h_cg
    - imr_.resulting_suspension_force_N.front_left * p_.b_f / 2
    + imr_.resulting_suspension_force_N.front_right * p_.b_f / 2
    - imr_.resulting_suspension_force_N.rear_left * p_.b_r / 2
    + imr_.resulting_suspension_force_N.rear_right * p_.b_r / 2
    + imr_.aero_torque_Nm.x
    + external_influences_.external_torque_Nm.x;
  if(p_.kinematics_enabled>.5) {
    const double h=p_.h_cg+x_vec_[x::z_m];
    imr_.resulting_torque_Nm.x=imr_.aero_torque_Nm.x+external_influences_.external_torque_Nm.x;
    imr_.resulting_torque_Nm.y=imr_.aero_torque_Nm.y-imr_.aero_force_N.z*(p_dep_.l_f_sprung-p_.l_f)+external_influences_.external_torque_Nm.y;
    imr_.resulting_torque_Nm.z=imr_.aero_torque_Nm.z+external_influences_.external_torque_Nm.z;
    for(int i=0;i<4;++i){const double side=i%2?-1.:1.,x=i<2?p_dep_.l_f_sprung:-p_dep_.l_r_sprung,y=side*(i<2?p_.b_f:p_.b_r)/2.-body_cg_y_m_,xplan=i<2?p_.l_f:-p_dep_.l_r;
      imr_.resulting_torque_Nm.x+=imr_.tire_forces_N[i][1]*h-imr_.resulting_suspension_force_N[i]*y;
      imr_.resulting_torque_Nm.y+=-imr_.tire_forces_N[i][0]*h+imr_.resulting_suspension_force_N[i]*x;
      imr_.resulting_torque_Nm.z+=imr_.tire_forces_N[i][1]*(xplan+kinematic_contact_x_[i])-imr_.tire_forces_N[i][0]*(kinematic_contact_y_[i]-vehicle_cg_y_m_);
    }
  }
  for(int i=0;i<4;++i){imr_.resulting_torque_Nm.x+=tire_mx_Nm_[i]*std::cos(imr_.effective_steering_angle_per_wheel_rad[i]);imr_.resulting_torque_Nm.y+=tire_mx_Nm_[i]*std::sin(imr_.effective_steering_angle_per_wheel_rad[i]);imr_.resulting_torque_Nm.z+=tire_mz_Nm_[i];}
  // clang-format on
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_drivetrain_load_torque_Nm()
{
  for (int i = 0; i < 4; ++i) {
    drivetrain_load_torque_per_wheel_Nm_[i] =
      (imr_.tire_forces_tire_frame_N[i][0] * contact_frames_[i].spin_lever_scale +
       imr_.tire_rolling_resistance_N[i]) *
      imr_.dynamic_tire_radius_m[i];
  }
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_steering_load_torque_Nm()
{
  // clang-format off
  #define CALC_SA_LOAD_TRQ(in_tire)\
    steering_load_torque_per_wheel_Nm_.in_tire = -tire_mz_Nm_.in_tire
    // Some fancy equation or model for torque steer
    // possibly from here: https://link.springer.com/content/pdf/10.1007/978-3-8348-9026-9_4

  EVAL_MACRO_PER_WHEEL(CALC_SA_LOAD_TRQ);
  #undef CALC_SA_LOAD_TRQ
  // clang-format on
}
}  // namespace tam::ocd::vehicle_dynamics
