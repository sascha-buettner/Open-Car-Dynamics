// Copyright 2026 Simon Sagmeister

#pragma once
#include <math.h>

#include <algorithm>
#include <array>
#include <eigen3/Eigen/Dense>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "ocd_aerodynamics_model_base_cpp/concept.hpp"
#include "ocd_tire_model_base_cpp/concept.hpp"
#include "ocd_types_cpp/road_surface.hpp"
#include "ocd_types_cpp/types.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/five_point_contact.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/mf62_vertical.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/nonlinear_radial_contact.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/road_contact_frame.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/states.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/suspension_lookup.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/suspension_table_2d.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/tire_load_memory.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/generic_radial_compliance.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/generic_radial_lookup.hpp"
#include "ocd_vehicle_dynamics_model_base_cpp/base_class.hpp"
#include "tum_helpers_cpp/constants.hpp"
#include "tum_helpers_cpp/numerical.hpp"
#include "tum_types_cpp/common.hpp"
namespace tam::ocd::vehicle_dynamics
{
// Forward declaration of the friend class
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
class VehicleDynamicsDoubleTrackModel;
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
class VehicleDynamicsDoubleTrackEqns
{
  friend class VehicleDynamicsDoubleTrackModel<TireModelT, AeroModelT>;
  // Types
  using x = tam::ocd::vehicle_dynamics::double_track::States::StateEnum;
  using StateVectorType = Eigen::Matrix<double, x::CNT_LENGTH_STATE_VECTOR, 1>;
  using double_per_wheel_t = tam::types::common::DataPerWheel<double>;
  using vector2d_per_wheel_t = tam::types::common::DataPerWheel<Eigen::Vector2d>;
  using vector3d_t = tam::types::common::Vector3D<double>;
  //
  struct Input
  {
    double_per_wheel_t wheel_speeds_radps{0, 0, 0, 0};       // Wheel speeds in rad/s
    double_per_wheel_t steering_angle_tire_rad{0, 0, 0, 0};  // Left is positive
  };
  struct Parameters
  {
    // Masses & Inertial moments
    double m;      // [kg]
    double m_t_f;  // [kg]
    double m_t_r;  // [kg]
    double I_z;    // [kg m^2]
    double I_y;    // [kg m^2]
    double I_x;    // [kg m^2]
    // Geometric parameters
    double l_f;     // [m]
    double l;       // [m]
    double b_f;     // [m]
    double b_r;     // [m]
    double h_cg;    // [m]
    double rr_w_f;  // [m]
    double rr_w_r;  // [m]
    // Height of the pitch center during acceleration
    double h_pc_accel;  // [m]
    // Height of the pitch center during acceleration
    double h_pc_decel;  // [m]
    // Height of the roll center
    tam::types::common::DataPerWheel<double> inflation_absolute_Pa{
      1.0, 1.0, 1.0, 1.0};
    double atmosphere_Pa{1.0};
    double h_rc;  // [m]
    double mf62_vertical_enabled{1.0}, mf62_rolling_enabled{1.0}, effective_radius_only{1.0};
    double longitudinal_relaxation_length_m{1.0};
    MF62Vertical vertical_front{}, vertical_rear{};
    double_per_wheel_t inflation_reference_absolute_Pa{1.0, 1.0, 1.0, 1.0};
    double kinematics_enabled{1.0}, corner_weights_enabled{1.0}, material_velocity_enabled{1.0};
    double_per_wheel_t corner_weights_kg{1.0, 1.0, 1.0, 1.0};
    std::vector<double> kinematic_travel{1.0, 1.0}, kinematic_mean_angle{1.0, 1.0},
      kinematic_rack{1.0, 1.0};
    SuspensionMap2D kinematic_front{}, kinematic_rear{};
    std::vector<double> rear_third_bar_travel{1.0, 1.0},
      rear_third_device_compression{1.0 / 1.0, 1.0 / 1.0};
    double geometry_enabled{1.0};
    double endstop_enabled{1.0}, endstop_stiffness{1.0}, endstop_damping{1.0};
    double stroke_bump_front{1.0}, stroke_rebound_front{1.0};
    double stroke_bump_rear{1.0}, stroke_rebound_rear{1.0};
    double caster_rad{1.0}, kingpin_inclination_rad{1.0};
    double h_rc_front{1.0}, h_rc_rear{1.0};
    double static_roll_centers_enabled{1.0};
    std::vector<double> geometry_travel{1.0, 1.0};
    std::vector<double> geometry_equiv_front{1.0, 1.0}, geometry_equiv_rear{1.0, 1.0};
    std::vector<double> geometry_camber_front{1.0, 1.0}, geometry_camber_rear{1.0, 1.0};
    std::vector<double> geometry_toe_front{1.0, 1.0}, geometry_toe_rear{1.0, 1.0};
    std::vector<double> geometry_roll_front{1.0, 1.0}, geometry_roll_rear{1.0, 1.0};
    std::vector<double> geometry_steer{1.0, 1.0, 1.0}, geometry_jacking{1.0, 1.0, 1.0};
    // Spring & Damping constants
    double c_f;  // [N/m]
    double c_r;  // [N/m]
    double d_f;  // [N*s/m]
    double d_r;  // [N*s/m]
    double damper_curve_enabled_f{1.0}, damper_curve_enabled_r{1.0};
    double damper_motion_ratio_f{1.0}, damper_motion_ratio_r{1.0};
    std::vector<double> damper_speed_mps{1.0, 1.0};
    std::vector<double> damper_bump_f{1.0, 1.0}, damper_rebound_f{1.0, 1.0};
    std::vector<double> damper_bump_r{1.0, 1.0}, damper_rebound_r{1.0, 1.0};
    // Optional radial-only Generic lookup. Tangential equations/states unchanged.
    double radial_lookup_enabled{1.0}, radial_conditioning_enabled{1.0}, radial_frozen_preload{1.0};
    double radial_stiffness_scale{1.0}, radial_nominal_stiffness{1.0};
    double low_load_boundary_N{1.0}, load_memory_fraction{1.0}, load_memory_tau_s{1.0};
    double fy_radial_mode{1.0}, fy_relaxation_strength{1.0};
    double radial_robust_solver_enabled{1.0};
    double force_conditioning_enabled{1.0};
    double_per_wheel_t force_condition_kappa{}, force_condition_alpha{}, radial_condition_fy_N{};
    double radial_speed_mps{1.0}, radial_softening_fraction{1.0};
    double radial_force_reference_f{1.0}, radial_force_reference_r{1.0};
    double_per_wheel_t radial_condition_load_N{1.0, 1.0, 1.0, 1.0};
    double_per_wheel_t radial_condition_fx_N{1.0, 1.0, 1.0, 1.0};
    double c_t_f;  // [N/m]
    double c_t_r;  // [N/m]
    // Anti-Roll bar
    double c_ar_f;  // [N/m]
    double c_ar_r;  // [N/m]
    // Local road contact: measured vertical/lateral compliance is configured per axle.
    double contact_enabled{1.0};
    double contact_sampling_points{1.0}, sidewall_enabled{1.0};
    double mf62_contact_force_coupling{1.0};
    double strict_surface_coverage{1.0};
    double tire_damping_f{1.0}, tire_damping_r{1.0};
    double resonant_belt_enabled{1.0}, resonant_belt_hz{1.0};
    double resonant_enabled{1.0}, resonant_hz{1.0}, resonant_width_hz{1.0}, resonant_residual{1.0};
    double maxwell_enabled{1.0}, maxwell_hz{1.0}, maxwell_fraction{1.0},
      maxwell_viscous_fraction{1.0};
    double_per_wheel_t resonant_acceleration{}, maxwell_derivative{};

    double tire_width_f{1.0}, tire_width_r{1.0};
    double contact_half_length{1.0};
    double contact_micro_radius{1.0};
    double enveloping_radius_f{1.0}, enveloping_radius_r{1.0};
    double radial_speed_override_enabled{1.0};
    double belt_enabled{1.0}, belt_mass_f{1.0}, belt_mass_r{1.0};
    double belt_frequency_f{1.0}, belt_frequency_r{1.0};
    double belt_damping_ratio_f{1.0}, belt_damping_ratio_r{1.0};
    double lateral_stiffness_f{1.0}, lateral_stiffness_r{1.0};
    double sidewall_height_f{1.0}, sidewall_height_r{1.0};
    double lateral_drop_scale{1.0};
    double lateral_drop_deflection_limit{1.0};
    double bump_gap_f{1.0}, bump_gap_r{1.0};
    double bump_stiffness_f{1.0}, bump_stiffness_r{1.0};
    // Optional third-element rubber. Device gap and wheel/device motion ratio.
    // Zero stiffness preserves the existing suspension; rubber curves require identification.
    double heave_stiffness_f{1.0}, heave_stiffness_r{1.0};
    double heave_spring_f{1.0}, heave_spring_r{1.0};
    double heave_gap_f{1.0}, heave_gap_r{1.0};
    double heave_motion_ratio_f{1.0}, heave_motion_ratio_r{1.0};
    // Rolling resistance
    double c_rr;  // https://en.wikipedia.org/wiki/Rolling_resistance
    // Tire Radii velocity scaling
    std::vector<double> rr_vel_scale_vel_points_mps;
    std::vector<double> rr_vel_scale_scale_factors;
    // Suspension geometry
    tam::types::common::DataPerWheel<double> toe_out_rad;
    tam::types::common::DataPerWheel<double> camber_out_rad;
  };
  // Variables that are directly calculated from the given parameters
  double_per_wheel_t radial_k_Npm_{1.0, 1.0, 1.0, 1.0};
  double_per_wheel_t radial_preload_N_{};
  double_per_wheel_t radial_force_coefficient_Npm_{};
  double_per_wheel_t radial_nominal_compression_m_{};
  std::array<GenericRadialCompliance, 4> radial_compliance_{};
  std::vector<double> radial_compliance_reference_signature_{};
  double radial_iteration_residual_{0.};
  double radial_root_fallback_used_{0.}, radial_root_evaluations_{0.};
  struct DependentParameters
  {
    double l_r = 0.0;
    double m_sprung = 0.0;    // Sprung mass
    double l_r_sprung = 0.0;  // Describing COG of the sprung mass
    double l_f_sprung = 0.0;  // Describing COG of the sprung mass
    // Coefficient for calculating anti squat and anti dive forces from the axle
    // See page A-30 of
    // https://babel.hathitrust.org/cgi/pt?id=mdp.39015075298698&view=1up&seq=146&skin=2021
    double slf = 0.0;
    double slr = 0.0;
    double sadf_deceleration = 0.0;
    double sadr_deceleration = 0.0;
    double sadf_acceleration = 0.0;
    double sadr_acceleration = 0.0;
  };
  struct IntermediateResults
  {
    // Sums up all effective steering angles including toe
    double_per_wheel_t effective_steering_angle_per_wheel_rad{0, 0, 0, 0};
    // Tire spring compression at standstill
    double_per_wheel_t tire_spring_initial_compression_m{0, 0, 0, 0};
    double_per_wheel_t tire_spring_force_N{0, 0, 0, 0};  // Tire spring forces
    double_per_wheel_t vertical_tire_force_N{0, 0, 0, 0};
    double_per_wheel_t antiroll_bar_force_N{0, 0, 0, 0};  // Antiroll bar forces
    // The actual velocity vector of the wheel in the contact patchoriginated in the vehicles
    // movement "Velocity over ground"
    vector2d_per_wheel_t velocity_wheel_over_ground_mps{
      Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero(),
      Eigen::Vector2d::Zero()};
    // The velocity vector of each tire rotated in each tires coordinate system
    vector2d_per_wheel_t velocity_wheel_over_ground_tire_frame_mps{
      Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero(),
      Eigen::Vector2d::Zero()};
    double_per_wheel_t velocity_tire_rotation_mps{0, 0, 0, 0};
    // Long slip and slip angle for each tire
    double_per_wheel_t tire_longitudinal_slip{0, 0, 0, 0};
    double_per_wheel_t tire_slip_angle_rad{0, 0, 0, 0};
    // Tire forces in longitudinal and lateral direction of each tire
    vector2d_per_wheel_t tire_forces_tire_frame_N{
      Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero(),
      Eigen::Vector2d::Zero()};
    vector2d_per_wheel_t tire_forces_N{
      Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero(),
      Eigen::Vector2d::Zero()};
    // Suspension spring compression at standstill
    double_per_wheel_t suspension_spring_initial_compression_m{0, 0, 0, 0};
    double_per_wheel_t suspension_spring_compression_m{0, 0, 0, 0};
    // Compression speed of the damper in mps
    double_per_wheel_t suspension_damper_compression_speed_mps{0, 0, 0, 0};
    // Forces in the suspension springs
    double_per_wheel_t suspension_spring_force_N{0, 0, 0, 0};
    // Forces in the dampers
    double_per_wheel_t suspension_damper_force_N{0, 0, 0, 0};
    // Torque on the tire because of rolling resistance
    // They won't be explicitly considered for calculating the lateral acceleration
    // since they will be considered as additional torque
    // keeping the tire from spinning faster
    double_per_wheel_t tire_rolling_resistance_N{0, 0, 0, 0};
    // Resulting vertical axle force
    double_per_wheel_t axle_vertical_force_N{0, 0, 0, 0};
    // Resulting force from suspension
    double_per_wheel_t resulting_suspension_force_N{0, 0, 0, 0};
    // Resulting veritcal force on the wheel
    double_per_wheel_t resulting_vertical_force_on_wheel_N{0, 0, 0, 0};
    // Aerodynamics
    vector3d_t aero_force_N{0, 0, 0};    // On the vehicle
    vector3d_t aero_torque_Nm{0, 0, 0};  // On the vehicle body
    // Resulting forces
    vector3d_t resulting_force_N{0, 0, 0};    // On the vehicle
    vector3d_t resulting_torque_Nm{0, 0, 0};  // On the vehicle body
    double_per_wheel_t road_contact_height_m{0, 0, 0, 0};
    // Diagnostic snapshots of the same five-node solve used by the dynamics.
    // Coordinates are in the moving planar road frame; z is a normal residual.
    std::array<double_per_wheel_t, 5> contact_node_query_x_m{}, contact_node_query_y_m{};
    std::array<double_per_wheel_t, 5> contact_node_road_z_m{}, contact_node_belt_z_m{};
    std::array<double_per_wheel_t, 5> contact_node_normal_N{};
    std::array<double_per_wheel_t, 5> contact_node_force_x_N{}, contact_node_force_y_N{};
    double_per_wheel_t contact_force_x_N{0, 0, 0, 0};
    double_per_wheel_t contact_force_y_N{0, 0, 0, 0};
    double_per_wheel_t contact_missing_samples{0, 0, 0, 0};
    double_per_wheel_t contact_bridged_samples{0, 0, 0, 0};
    double_per_wheel_t tire_lateral_deflection_m{0, 0, 0, 0};
    double_per_wheel_t tire_lateral_drop_m{0, 0, 0, 0};
    double_per_wheel_t tire_lateral_target_m{0, 0, 0, 0};
    double_per_wheel_t tire_deformation_extrapolated{0, 0, 0, 0};
    double_per_wheel_t suspension_bump_force_N{0, 0, 0, 0};
    double_per_wheel_t endstop_force_N{0, 0, 0, 0};
    double_per_wheel_t damper_curve_speed_extrapolated{0, 0, 0, 0};
    double_per_wheel_t suspension_heave_force_N{0, 0, 0, 0};
    double_per_wheel_t suspension_heave_energy_J{0, 0, 0, 0};
    double_per_wheel_t contact_yaw_moment_Nm{0, 0, 0, 0};
    double_per_wheel_t contact_elastic_energy_J{0, 0, 0, 0};
    double_per_wheel_t belt_carcass_force_N{0, 0, 0, 0};
    double_per_wheel_t belt_resulting_force_N{0, 0, 0, 0};
    double_per_wheel_t belt_elastic_energy_J{0, 0, 0, 0};
    double_per_wheel_t rim_support_force_N{0, 0, 0, 0};
    double_per_wheel_t contact_normal_load_N{0, 0, 0, 0};
    double_per_wheel_t contact_tangential_vertical_force_N{0, 0, 0, 0};
    // Dynamic Tire Radius
    double_per_wheel_t dynamic_tire_radius_m{0, 0, 0, 0};
  };
  //
  std::vector<double> radial_reference_signature_;
  double_per_wheel_t radial_reference_loaded_m_{0, 0, 0, 0};
  double_per_wheel_t mf62_free_radius_m_{0, 0, 0, 0}, mf62_loaded_radius_m_{0, 0, 0, 0},
    mf62_deflection_m_{0, 0, 0, 0};
  double_per_wheel_t kinematic_arb_m_{0, 0, 0, 0}, kinematic_arb_jac_{1, 1, 1, 1};
  double_per_wheel_t tire_mx_Nm_{0, 0, 0, 0}, tire_mz_Nm_{0, 0, 0, 0};
  double_per_wheel_t tire_gamma_rad_{0, 0, 0, 0}, radial_residual_N_{0, 0, 0, 0};
  double body_cg_y_m_{0.}, vehicle_cg_y_m_{0.}, kinematic_rack_m_{0.}, kinematic_rack_rate_mps_{0.};
  double_per_wheel_t kinematic_q_m_{0, 0, 0, 0}, kinematic_q_velocity_mps_{0, 0, 0, 0};
  double_per_wheel_t kinematic_spring_m_{0, 0, 0, 0}, kinematic_spring_jac_{1, 1, 1, 1},
    kinematic_spring_steer_jac_{0, 0, 0, 0}, kinematic_damper_steer_jac_{0, 0, 0, 0};
  double_per_wheel_t kinematic_jac_x_{0, 0, 0, 0}, kinematic_jac_y_{0, 0, 0, 0},
    kinematic_jac_z_{1, 1, 1, 1};
  double_per_wheel_t kinematic_contact_x_{0, 0, 0, 0}, kinematic_contact_y_{0, 0, 0, 0},
    kinematic_contact_z_{0, 0, 0, 0};
  double_per_wheel_t kinematic_angjac_x_{0, 0, 0, 0}, kinematic_angjac_y_{0, 0, 0, 0},
    kinematic_angjac_z_{0, 0, 0, 0};
  double roll_center_front_y_m_{0.}, roll_center_rear_y_m_{0.};
  double steering_mean_rate_radps_{0.};
  double_per_wheel_t geometry_equivalent_extension_m_{0, 0, 0, 0}, geometry_leverage_{1, 1, 1, 1};
  double_per_wheel_t kinematic_jacking_m_{0, 0, 0, 0};
  double_per_wheel_t kinematic_camber_rad_{0, 0, 0, 0};
  double_per_wheel_t kinematic_toe_rad_{0, 0, 0, 0};
  double_per_wheel_t wheel_speeds_radps_{0, 0, 0, 0};
  double_per_wheel_t steering_angle_per_wheel_rad_{0, 0, 0, 0};
  //
  StateVectorType x_vec_ = StateVectorType::Zero();
  StateVectorType x_dot_vec_ = StateVectorType::Zero();
  //
  double_per_wheel_t drivetrain_load_torque_per_wheel_Nm_{0, 0, 0, 0};
  double_per_wheel_t steering_load_torque_per_wheel_Nm_{0, 0, 0, 0};
  types::VehicleDynamicsModelOutput vd_output_{};
  //
  IntermediateResults imr_{};
  std::array<RoadContactFrame, 4> contact_frames_{};
  //
  types::ExternalInfluences external_influences_{};
  std::shared_ptr<const types::RoadSurface> road_surface_;
  double_per_wheel_t point_road_gx_{0, 0, 0, 0}, point_road_gy_{0, 0, 0, 0},
    point_road_h_{0, 0, 0, 0};
  FivePointContact contact_;
  void calculate_local_road_contact();
  tam::types::common::DataPerWheel<TireModelT> tire_models_;
  AeroModelT aero_model_;
  Parameters p_;
  DependentParameters p_dep_{};
  //
  // Constants
  static constexpr double G_CONST = 9.81;  // mps2
  double_per_wheel_t tire_spring_initial_compression_{0, 0, 0, 0};
  //
  double radial_lookup_speed_mps() const
  {
    return p_.radial_conditioning_enabled > .5 || p_.radial_speed_override_enabled > .5
             ? p_.radial_speed_mps
             : std::hypot(x_vec_[x::v_x_mps], x_vec_[x::v_y_mps]);
  }
  void calculate_intermediate_results();
  void calculate_suspension_output();
  // Function to calculate pars of the intermediate results.
  void calc_dynamic_tire_radius();
  void calculate_effective_steering_angle();
  void calculate_dependent_parameters();
  void calculate_tire_spring_initial_compression_m();
  void calculate_tire_spring_force_N();
  void calculate_vertical_tire_force_N();
  void calculate_velocity_wheel_over_ground_mps();
  void calculate_velocity_wheel_over_ground_tire_frame_mps();
  void calculate_velocity_tire_rotation_mps();
  void calculate_tire_longitudinal_slip();
  void calculate_tire_slip_angle_rad();
  void calculate_tire_forces_tire_frame_N();
  void calculate_tire_forces_N();
  void calculate_suspension_spring_initial_compression_m();
  void calculate_suspension_spring_compression_m();
  void calculate_suspension_damper_compression_speed_mps();
  void calculate_suspension_spring_force_N();
  void calculate_suspension_damper_force_N();
  void calculate_antiroll_bar_force_N();
  void calculate_tire_rolling_resistance_N();
  void calculate_axle_vertical_force_N();
  void calculate_resulting_suspension_force_N();
  void calculate_resulting_vertical_force_on_wheel_N();
  void calculate_aerodynamics();
  void calculate_resulting_force_N();
  void calculate_resulting_torque_Nm();
  void calculate_drivetrain_load_torque_Nm();
  void calculate_steering_load_torque_Nm();
  //
public:
  void set_wheel_speeds(const double_per_wheel_t & wheel_speeds_radps)
  {
    wheel_speeds_radps_ = wheel_speeds_radps;
  }
  void set_steering_rate(double rate) { steering_mean_rate_radps_ = rate; }
  void set_steering_angles(const double_per_wheel_t & steering_angle_per_wheel_rad)
  {
    steering_angle_per_wheel_rad_ = steering_angle_per_wheel_rad;
  }
  void set_external_influences(const types::ExternalInfluences & external_influences)
  {
    external_influences_ = external_influences;
  }
  void set_x_vec(const StateVectorType & x_vec) { x_vec_ = x_vec; }
  //
  // Evaluate the ODE
  void evaluate();
};
}  // namespace tam::ocd::vehicle_dynamics
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_equations_impl.hpp"
