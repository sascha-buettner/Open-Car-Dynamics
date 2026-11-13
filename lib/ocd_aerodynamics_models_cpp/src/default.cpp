// Copyright 2026 Simon Sagmeister
#include "ocd_aerodynamics_models_cpp/default.hpp"

#include <cmath>

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "tum_helpers_cpp/numerical.hpp"
namespace tam::ocd::aerodynamics
{
types::AeroModelOutput DefaultAerodynamicsModel::evaluate(types::AeroModelInput const & u)
{
  // vx, vy - long/lat velocity in m/s, theta - pitch angle in rad
  auto sign = [](double x) {
    if (x < 0) {
      return -1;
    } else {
      return 1;
    }
  };

  types::AeroModelOutput out;
  double const vx_rel_mps = u.vx_mps - u.wind_mps.x;
  double const vy_rel_mps = u.vy_mps - u.wind_mps.y;

  // Lift cofficient multiplier interpolation
  double delta_c_l = tam::helpers::numerical::interp(
    vx_rel_mps, p_.c_l_velocity_scaling__vel_mps, p_.c_l_velocity_scaling__delta_c_l);
  // Aero center translation interpolation
  double aero_center_translation_m = tam::helpers::numerical::interp(
    u.pitch_angle_rad, p_.pitch_aero_center_translation__pitch_rad,
    p_.pitch_aero_center_translation__aero_center_translation_m);

  out.force_cog_N.x =
    -0.5 * p_.c_d * p_.air_density_kgpm3 * p_.A_m2 * std::pow(vx_rel_mps, 2) * sign(vx_rel_mps);
  // To be precise the area of calculating the drag here is different for lateral movement
  // but this is neglected here
  out.force_cog_N.y =
    -0.5 * p_.c_d * p_.air_density_kgpm3 * p_.A_m2 * std::pow(vy_rel_mps, 2) * sign(vy_rel_mps);
  out.force_cog_N.z =
    0.5 * (p_.c_l + delta_c_l) * p_.air_density_kgpm3 * p_.A_m2 * std::pow(vx_rel_mps, 2);

  if (p_.attitude_enabled > .5) {
    if (!(std::isfinite(p_.attitude_wheelbase_m) && p_.attitude_wheelbase_m > 0. &&
          p_.attitude_front_distance_m > 0. &&
          p_.attitude_front_distance_m < p_.attitude_wheelbase_m &&
          p_.attitude_minimum_speed_mps > 0.))
      throw std::invalid_argument("Invalid attitude aero geometry");
    auto bound = [](double x, double lo, double hi) {
      if (!(lo < 0. && hi > 0.)) throw std::invalid_argument("Invalid attitude aero bounds");
      return x >= 0. ? hi * std::tanh(x / hi) : -lo * std::tanh(x / (-lo));
    };
    // Our theta is nose-down positive. Marchesin et al. Eq21 uses the
    // opposite pitch sign: theta_eff=theta+L*q/U-vz/U in our coordinates.
    const double U = std::max(std::abs(vx_rel_mps), p_.attitude_minimum_speed_mps);
    const double theta =
      u.pitch_angle_rad + p_.attitude_unsteady_gain *
                            (p_.attitude_wheelbase_m * u.pitch_rate_radps - u.heave_velocity_mps) /
                            U;
    const double dz = u.z_m - p_.attitude_reference_z_m;
    const double hf = dz - p_.attitude_front_distance_m * theta;
    const double hr = dz + (p_.attitude_wheelbase_m - p_.attitude_front_distance_m) * theta;
    const double lift = bound(
      p_.attitude_lift_front_gradient_per_m * hf + p_.attitude_lift_rear_gradient_per_m * hr,
      p_.attitude_lift_delta_min, p_.attitude_lift_delta_max);
    const double balance = bound(
      p_.attitude_balance_front_gradient_per_m * hf + p_.attitude_balance_rear_gradient_per_m * hr,
      p_.attitude_balance_delta_min, p_.attitude_balance_delta_max);
    out.force_cog_N.z *= 1. + p_.attitude_lift_gain * lift;
    aero_center_translation_m += p_.attitude_center_gain * p_.attitude_wheelbase_m * balance;
  }
  // Aero Balance
  out.torque_Nm.y = -out.force_cog_N.z * (p_.d_aero_center_cog_m + aero_center_translation_m);

  // Zero for everything else
  out.torque_Nm.x = 0;
  out.torque_Nm.z = 0;

  return out;
}
void DefaultAerodynamicsModel::declare_parameters(
  tam::pmg::ParamReferenceManager * param_manager, std::string name_prefix)
{
  param_manager->declare_parameter(
    name_prefix + "air_density_kgpm3", &(p_.air_density_kgpm3), 1.225,
    tam::pmg::ParameterType::DOUBLE, "");

  // Base parameters
  param_manager->declare_parameter(
    name_prefix + "A_m2", &(p_.A_m2), 1.0, tam::pmg::ParameterType::DOUBLE, "");
  param_manager->declare_parameter(
    name_prefix + "c_d", &(p_.c_d), 1.0, tam::pmg::ParameterType::DOUBLE, "");
  param_manager->declare_parameter(
    name_prefix + "c_l", &(p_.c_l), -1.5, tam::pmg::ParameterType::DOUBLE, "");
  param_manager->declare_parameter(
    name_prefix + "d_aero_center_cog_m", &(p_.d_aero_center_cog_m), 0.0,
    tam::pmg::ParameterType::DOUBLE, "");

  param_manager->declare_parameter(
    name_prefix + "attitude_map.enabled", &p_.attitude_enabled, 1.0,
    tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.lift_gain", &p_.attitude_lift_gain, 1.0,
    tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.center_gain", &p_.attitude_center_gain, 1.0,
    tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.unsteady_gain", &p_.attitude_unsteady_gain, 1.0,
    tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.reference_z_m", &p_.attitude_reference_z_m, 1.0,
    tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.wheelbase_m", &p_.attitude_wheelbase_m, 1.0,
    tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.front_distance_m", &p_.attitude_front_distance_m, 1.0, tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.lift_front_gradient_per_m", &p_.attitude_lift_front_gradient_per_m, 1.0, tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.lift_rear_gradient_per_m", &p_.attitude_lift_rear_gradient_per_m, 1.0, tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.balance_front_gradient_per_m",
    &p_.attitude_balance_front_gradient_per_m, 1.0, tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.balance_rear_gradient_per_m",
    &p_.attitude_balance_rear_gradient_per_m, 1.0, tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.lift_delta_min", &p_.attitude_lift_delta_min, 1.0,
    tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.lift_delta_max", &p_.attitude_lift_delta_max, 1.0,
    tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.balance_delta_min", &p_.attitude_balance_delta_min, 1.0,
    tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.balance_delta_max", &p_.attitude_balance_delta_max, 1.0,
    tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");
  param_manager->declare_parameter(
    name_prefix + "attitude_map.minimum_speed_mps", &p_.attitude_minimum_speed_mps, 1.0,
    tam::pmg::ParameterType::DOUBLE,
    "Bounded literature-normalized sensitivity; not an identified DummyVehicle map");

  // Velocity scaling of cl
  param_manager->declare_parameter(
    name_prefix + "c_l_velocity_scaling.vel_mps", &(p_.c_l_velocity_scaling__vel_mps),
    std::vector<double>{0, 100}, tam::pmg::ParameterType::DOUBLE_ARRAY, "");

  param_manager->declare_parameter(
    name_prefix + "c_l_velocity_scaling.delta_c_l", &(p_.c_l_velocity_scaling__delta_c_l),
    std::vector<double>{0, 0}, tam::pmg::ParameterType::DOUBLE_ARRAY, "");

  // Pitch correction since pitch moves the center of aero forces
  param_manager->declare_parameter(
    name_prefix + "pitch_aero_center_translation.pitch_rad",
    &(p_.pitch_aero_center_translation__pitch_rad), std::vector<double>{-M_PI, M_PI},
    tam::pmg::ParameterType::DOUBLE_ARRAY, "");

  param_manager->declare_parameter(
    name_prefix + "pitch_aero_center_translation.aero_center_translation_m",
    &(p_.pitch_aero_center_translation__aero_center_translation_m), std::vector<double>{0, 0},
    tam::pmg::ParameterType::DOUBLE_ARRAY, "");
}
}  // namespace tam::ocd::aerodynamics
