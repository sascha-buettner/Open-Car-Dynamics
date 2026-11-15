// Copyright 2026 TUM Autonomous Motorsport
#pragma once
#include <eigen3/Eigen/Cholesky>

namespace tam::ocd::vehicle_dynamics {
// Open, small-displacement six-coordinate ring, in nonrotating hub axes.
// Absolute belt velocities avoid differentiating delayed spindle inputs.
// Vertical translation uses the existing mass-partitioned belt/contact pair.
template <tam::ocd::interfaces::concepts::TireModel TireModelT,
          tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT,AeroModelT>::calculate_ring_reactions()
{
  ring_pitch_correction_.setZero();
  if (p_.ring_enabled < .5) return;
  if (p_.belt_enabled < .5 || p_.contact_enabled < .5)
    throw std::invalid_argument("mfSCHWIFT needs explicit ring radial/contact dynamics");
  for(double value:{p_.ring_inertia_trans,p_.ring_inertia_spin,p_.ring_frequency_lat,p_.ring_frequency_yaw,p_.ring_frequency_windup,p_.ring_longitudinal_stiffness})
    if(!(std::isfinite(value)&&value>0.))throw std::invalid_argument("mfSCHWIFT positive structural parameter required");
  for(double value:{p_.ring_damping_lat,p_.ring_damping_yaw,p_.ring_damping_windup})
    if(!(std::isfinite(value)&&value>=0.))throw std::invalid_argument("mfSCHWIFT nonnegative damping required");
  if(!(p_.ring_lateral_shear_fraction>=0.&&p_.ring_lateral_shear_fraction<=1.))throw std::invalid_argument("mfSCHWIFT shear centre outside declared radius bounds");
  if(p_.maxwell_enabled>.5||p_.resonant_belt_enabled>.5)throw std::invalid_argument("mfSCHWIFT modal sidewall cannot silently inherit empirical belt rheology");
  for (int i=0;i<4;++i) {
    const double psi=x_vec_[x::ring_yaw_FL_rad+i];
    const double delta=imr_.effective_steering_angle_per_wheel_rad[i]-psi;
    const double c=std::cos(delta),s=std::sin(delta),mb=i<2?p_.belt_mass_f:p_.belt_mass_r;
    const double friction_x=imr_.tire_forces_N[i][0],friction_y=imr_.tire_forces_N[i][1];
    const double fx=friction_x+imr_.contact_force_x_N[i],fy=friction_y+imr_.contact_force_y_N[i];
    const double Fxt=c*friction_x+s*friction_y,Fyt=-s*friction_x+c*friction_y,R=imr_.dynamic_tire_radius_m[i];
    const double qx=x_vec_[x::ring_x_FL_m+i],qy=x_vec_[x::ring_y_FL_m+i];
    const double xv=x_vec_[x::ring_vx_FL_mps+i],yv=x_vec_[x::ring_vy_FL_mps+i];
    const double wx=i<2?p_.l_f:-p_dep_.l_r,
                 wy=(i%2?-1.:1.)*(i<2?p_.b_f:p_.b_r)/2.;
    const double r=x_vec_[x::psi_dot_radps];
    const double joint_rate=p_.kinematics_enabled>.5?kinematic_q_velocity_mps_[i]:0.;
    const double hub_p=x_vec_[x::phi_dot_radps]+(p_.kinematics_enabled>.5?kinematic_angjac_x_[i]*joint_rate:0.);
    const double hub_q=x_vec_[x::theta_dot_radps]+(p_.kinematics_enabled>.5?kinematic_angjac_y_[i]*joint_rate:0.);
    const double hub_r=r+(p_.kinematics_enabled>.5?kinematic_angjac_z_[i]*joint_rate:0.)+(i<2?steering_mean_rate_radps_:0.);
    const double hub_vx=imr_.velocity_wheel_over_ground_mps[i][0]+R*hub_q;
    const double hub_vy=imr_.velocity_wheel_over_ground_mps[i][1]-R*hub_p;
    const double relx=xv-hub_vx+r*qy;
    const double rely=yv-hub_vy-r*qx;
    const double f_long=i<2?p_.belt_frequency_f:p_.belt_frequency_r;
    const double kz=mb*std::pow(2*M_PI*f_long,2),ky=mb*std::pow(2*M_PI*p_.ring_frequency_lat,2);
    const double dx=2*(i<2?p_.belt_damping_ratio_f:p_.belt_damping_ratio_r)*mb*2*M_PI*f_long;
    const double dy=2*p_.ring_damping_lat*mb*2*M_PI*p_.ring_frequency_lat;
    const double Sxt=kz*(c*qx+s*qy)+dx*(c*relx+s*rely);
    const double shear_lever=p_.ring_lateral_shear_fraction*R;
    const double shear_gamma=x_vec_[x::ring_camber_FL_rad+i];
    const double shear_rate=x_vec_[x::ring_omega_x_FL_radps+i]-(c*hub_p+s*hub_q);
    const double Syt=ky*(-s*qx+c*qy+shear_lever*shear_gamma)+dy*(-s*relx+c*rely+shear_lever*shear_rate);
    const double Sx=c*Sxt-s*Syt,Sy=s*Sxt+c*Syt;
    ring_derivatives_[0][i]=relx;ring_derivatives_[1][i]=rely;
    ring_derivatives_[2][i]=(fx-Sx)/mb+r*yv;
    ring_derivatives_[3][i]=(fy-Sy)/mb-r*xv;
    const double qg=x_vec_[x::ring_camber_FL_rad+i],qt=x_vec_[x::ring_windup_FL_rad+i],
                 qp=x_vec_[x::ring_yaw_FL_rad+i];
    const double bg=x_vec_[x::ring_omega_x_FL_radps+i],bt=x_vec_[x::ring_omega_y_FL_radps+i],
                 bp=x_vec_[x::ring_omega_z_FL_radps+i];
    const double hg=c*hub_p+s*hub_q;
    const double ht=-s*hub_p+c*hub_q+wheel_speeds_radps_[i];
    const double hp=hub_r;
    const double vg=bg-hg,vt=bt-ht,vp=bp-hp;
    const double kg=p_.ring_inertia_trans*std::pow(2*M_PI*p_.ring_frequency_yaw,2);
    const double kt=p_.ring_inertia_spin*std::pow(2*M_PI*p_.ring_frequency_windup,2);
    const double dg=2*p_.ring_damping_yaw*p_.ring_inertia_trans*2*M_PI*p_.ring_frequency_yaw;
    const double dt=2*p_.ring_damping_windup*p_.ring_inertia_spin*2*M_PI*p_.ring_frequency_windup;
    const double Sg=kg*qg+dg*vg+shear_lever*Syt,St=kt*qt+dt*vt,Sp=kg*qp+dg*vp;
    double normal_mx=R*(-s*imr_.contact_force_x_N[i]+c*imr_.contact_force_y_N[i]);
    if (p_.contact_sampling_points==5.) {
      const double width=(i<2?p_.tire_width_f:p_.tire_width_r)/2.;
            const double py[5]={0.,width,-width,width,-width};
      for(int j=0;j<5;++j){normal_mx+=py[j]*imr_.contact_node_normal_N[j][i];}
    }
    const double Mg=Fyt*R+tire_mx_Nm_[i]+normal_mx;
    // A circular ring is invariant under spin: the conservative normal
    // envelope has no windup-coordinate torque. Tangential traction and MF
    // rolling resistance act on spin; normal road gradients act on translation.
    const double Mt=-Fxt*R-imr_.tire_rolling_resistance_N[i]*R*std::tanh(bt*R/.1);
    const double Mp=tire_mz_Nm_[i]+imr_.contact_yaw_moment_Nm[i];
    // Gyroscopic matrix is skew: it does no work. This open modal-support
    // approximation does not reproduce proprietary rotating-damper corrections.
    const double gyro=p_.ring_inertia_spin*bt;
    ring_derivatives_[4][i]=vg;ring_derivatives_[5][i]=vt;ring_derivatives_[6][i]=vp;
    ring_derivatives_[7][i]=(Mg-Sg+gyro*bp)/p_.ring_inertia_trans;
    ring_derivatives_[8][i]=(Mt-St)/p_.ring_inertia_spin;
    ring_derivatives_[9][i]=(Mp-Sp-gyro*bg)/p_.ring_inertia_trans;
    ring_spin_load_Nm_[i]=-St;
    const double nominal=p_.radial_lookup_enabled>.5?radial_preload_N_[i]:imr_.tire_spring_initial_compression_m[i]*(i<2?p_.c_t_f:p_.c_t_r);
    ring_radial_spindle_N_[i]=kz*(x_vec_[x::z_belt_FL_m+i]-x_vec_[x::z_w_FL_m+i]+(nominal-mb*G_CONST)/kz)
       +dx*(x_vec_[x::v_z_belt_FL_mps+i]-x_vec_[x::v_z_w_FL_mps+i]);
    const double gamma_correction=Sg-Syt*R,pitch_correction=St+Sxt*R;
    ring_tire_pitch_correction_Nm_[i]=-std::sin(psi)*gamma_correction+std::cos(psi)*pitch_correction;
    // Spindle force is sidewall reaction, never a second contact force.
    imr_.tire_forces_N[i]={Sx,Sy};
    tire_mx_Nm_[i]=std::cos(psi)*gamma_correction+std::sin(psi)*pitch_correction;
    tire_mz_Nm_[i]=Sp+qx*Sy-qy*Sx;
    imr_.contact_yaw_moment_Nm[i]=0.;
    ring_pitch_correction_[0]-=std::sin(delta+psi)*ring_tire_pitch_correction_Nm_[i];
    ring_pitch_correction_[1]+=std::cos(delta+psi)*ring_tire_pitch_correction_Nm_[i];
    ring_stored_energy_[i]=.5*mb*(xv*xv+yv*yv)+.5*kz*std::pow(c*qx+s*qy,2)+.5*ky*std::pow(-s*qx+c*qy+shear_lever*qg,2)
       +.5*p_.ring_inertia_trans*(bg*bg+bp*bp)+.5*p_.ring_inertia_spin*bt*bt
       +.5*kg*(qg*qg+qp*qp)+.5*kt*qt*qt;
    ring_loss_[i]=dx*std::pow(c*relx+s*rely,2)+dy*std::pow(-s*relx+c*rely+shear_lever*vg,2)+dg*(vg*vg+vp*vp)+dt*vt*vt;
  }
}

template <tam::ocd::interfaces::concepts::TireModel TireModelT,
          tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
Eigen::Vector3d VehicleDynamicsDoubleTrackEqns<TireModelT,AeroModelT>::ring_planar_acceleration() const
{
  Eigen::Matrix3d mass=Eigen::Vector3d(p_.m,p_.m,p_.I_z).asDiagonal();
  Eigen::Vector3d rhs(imr_.resulting_force_N.x,imr_.resulting_force_N.y,imr_.resulting_torque_Nm.z);
  Eigen::Vector2d first=Eigen::Vector2d::Zero();
  for(int i=0;i<4;++i){
    const double mb=i<2?p_.belt_mass_f:p_.belt_mass_r,xx=i<2?p_.l_f:-p_dep_.l_r,
                 yy=(i%2?-1.:1.)*(i<2?p_.b_f:p_.b_r)/2.;
    Eigen::Matrix<double,2,3> J;J<<1.,0.,-yy,0.,1.,xx;
    mass.noalias()-=mb*J.transpose()*J;mass(2,2)-=p_.ring_inertia_trans;
    first-=mb*Eigen::Vector2d(xx,yy);
  }
  rhs.head<2>()+=std::pow(x_vec_[x::psi_dot_radps],2)*first;
  Eigen::LDLT<Eigen::Matrix3d> solver(mass);
  if(!solver.isPositive())throw std::domain_error("mfSCHWIFT mass partition is not positive definite");
  return solver.solve(rhs);
}
} // namespace tam::ocd::vehicle_dynamics
