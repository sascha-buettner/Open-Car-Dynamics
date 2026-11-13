// Copyright 2026 Simon Sagmeister
#pragma once
#include <array>
#include <vector>

#include "tum_types_cpp/common.hpp"
namespace tam::ocd::types
{
struct ExternalInfluences
{
  tam::types::common::Vector3D<double> external_force_N = {0, 0, 0};    // On the vehicle body
  tam::types::common::Vector3D<double> external_torque_Nm = {0, 0, 0};  // On the vehicle body
  tam::types::common::Vector3D<double> wind_mps = {0, 0, 0};            // In vehicle frame
  tam::types::common::DataPerWheel<double> z_height_road_m = {0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> lambda_mue = {1, 1, 1, 1};
  // Diagnostic vertical forces, positive upwards in the yaw-aligned road frame.
  // Rim force acts on the unsprung rim mass. Tire force acts on the belt mass
  // when enabled, otherwise on the combined unsprung mass. No body wrench is
  // added a second time: these loads reach the body through the suspension.
  tam::types::common::DataPerWheel<double> rim_force_z_N = {0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> tire_force_z_N = {0, 0, 0, 0};
};
// Read-only live diagnostic snapshot. This is not a cached logger data map and
// is intentionally separate from the transport/public observation layout.
struct TireContactDiagnostics
{
  tam::types::common::DataPerWheel<double> tire_spring_force_N{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> belt_displacement_m{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> belt_velocity_mps{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> dynamic_radius_m{0, 0, 0, 0};
  std::array<tam::types::common::DataPerWheel<double>, 5> node_normal_N{};
};
// Generalized sprung-body motion relative to the yaw-aligned road plane.
// Angles are Euler roll/pitch; their derivatives are not body-frame gyro rates.
struct BodyMotion
{
  double heave_m{0.0}, heave_velocity_mps{0.0}, heave_acceleration_mps2{0.0};
  double roll_rad{0.0}, roll_rate_radps{0.0}, roll_acceleration_radps2{0.0};
  double pitch_rad{0.0}, pitch_rate_radps{0.0}, pitch_acceleration_radps2{0.0};
};
struct SuspensionOutput
{
  tam::types::common::DataPerWheel<double> third_equivalent_extension_m{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> damper_extension_m{0, 0, 0, 0};
  // Vertical displacements from the model's nominal geometry, positive upwards.
  // Road-relative values subtract ExternalInfluences::z_height_road_m per wheel.
  tam::types::common::DataPerWheel<double> body_displacement_m{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> wheel_displacement_m{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> body_displacement_over_road_m{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> wheel_displacement_over_road_m{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> wheel_velocity_mps{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> wheel_acceleration_mps2{0, 0, 0, 0};
  // Existing model spring/damper convention: compression and compressive forces are negative.
  tam::types::common::DataPerWheel<double> spring_compression_m{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> damper_compression_velocity_mps{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> spring_force_N{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> damper_force_N{0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> bump_force_N{0, 0, 0, 0};
};
struct VehicleDynamicsModelOutput
{
  // Planar motion supplied to the road-geometry coupling.
  tam::types::common::Vector3D<double> position_m = {0, 0, 0};
  tam::types::common::Vector3D<double> velocity_mps = {0, 0, 0};
  tam::types::common::Vector3D<double> acceleration_mps2 = {0, 0, 0};
  tam::types::common::Vector3D<double> orientation_rad = {0, 0, 0};
  tam::types::common::Vector3D<double> angular_velocity_radps = {0, 0, 0};
  tam::types::common::Vector3D<double> angular_acceleration_radps2 = {0, 0, 0};
  tam::types::common::DataPerWheel<double> tire_longitudinal_slip_per_wheel = {0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> tire_slip_angle_per_wheel_rad = {0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> longitudinal_tire_force_tire_frame_per_wheel_N = {
    0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> lateral_tire_force_tire_frame_per_wheel_N = {0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> vertical_tire_force_per_wheel_N = {0, 0, 0, 0};
  BodyMotion body_motion{};
  SuspensionOutput suspension{};
};
struct VehicleModelOutput
{
  VehicleDynamicsModelOutput vehicle_dynamics_output;
  tam::types::common::DataPerWheel<double> wheel_speeds_radps = {0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> drivetrain_load_torque_per_wheel_Nm = {0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> steering_angle_per_wheel_rad = {0, 0, 0, 0};
  tam::types::common::DataPerWheel<double> steering_load_torque_per_wheel_Nm = {0, 0, 0, 0};
};
struct AeroModelOutput
{
  tam::types::common::Vector3D<double> force_cog_N = {0, 0, 0};
  tam::types::common::Vector3D<double> torque_Nm = {0, 0, 0};
};
struct AeroModelInput
{
  double vx_mps = 0.0;                                        // Velocity in x direction
  double vy_mps = 0.0;                                        // Velocity in y direction
  tam::types::common::Vector3D<double> wind_mps = {0, 0, 0};  // Air velocity in vehicle frame
  double pitch_angle_rad = 0.0;                               // Pitch angle of the vehicle
  double z_m = 0.0;                                           // Heave from neutral position
  double pitch_rate_radps = 0.0;
  double heave_velocity_mps = 0.0;
};
struct TireModelOutput
{
  double longitudinal_force_N = 0.0;  // Tire force in longitudinal direction
  double lateral_force_N = 0.0;       // Tire force in lateral direction
  double overturning_moment_Nm = 0.0;
  double self_aligning_moment_Nm = 0.0;  // Self aligning moment
};
struct TireModelInput
{
  double longitudinal_slip = 0.0;  // Slip ratio in longitudinal direction
  double slip_angle_rad = 0.0;     // Slip angle at the respective tire
  double F_z_N = 0.0;              // Normal force in the tire road contact patch
  double gamma_rad =
    0.0;  // Positive angle means clockwise rotation locking from the rear of the car.
  double inflation_pressure_Pa = 1.0;  // Inert placeholder; configure explicitly
};
}  // namespace tam::ocd::types
