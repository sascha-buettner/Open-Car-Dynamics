#include <cmath>
#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
using tam::ocd::aerodynamics::DefaultAerodynamicsModel;

int main()
{
  DefaultAerodynamicsModel model;
  auto & p = model.p_;
  p.air_density_kgpm3 = 1.2;
  p.c_d = .8;
  p.c_l = 1.5;
  p.A_m2 = 1.4;
  p.d_aero_center_cog_m = .1;
  p.c_l_velocity_scaling__vel_mps = {0., 50.};
  p.c_l_velocity_scaling__delta_c_l = {0., 0.};
  p.pitch_aero_center_translation__pitch_rad = {-.2, 0., .2};
  p.pitch_aero_center_translation__aero_center_translation_m = {0., 0., 0.};
  tam::ocd::types::AeroModelInput input;
  input.vx_mps = 30.;
  const auto plain = model.evaluate(input);
  const double force = .5 * p.c_l * p.air_density_kgpm3 * p.A_m2 * 30. * 30.;
  if (std::abs(plain.force_cog_N.z - force) > 1e-9 ||
      std::abs(plain.torque_Nm.y + force * p.d_aero_center_cog_m) > 1e-9)
    throw std::runtime_error("Public fixed aero baseline changed");

  p.attitude_enabled = 1.;
  p.attitude_lift_gain = p.attitude_center_gain = 1.;
  p.attitude_unsteady_gain = 1.;
  p.attitude_reference_z_m = 0.;
  p.attitude_wheelbase_m = 3.;
  p.attitude_front_distance_m = 1.4;
  p.attitude_minimum_speed_mps = 5.;
  p.attitude_lift_front_gradient_per_m = .2;
  p.attitude_lift_rear_gradient_per_m = .05;
  p.attitude_balance_front_gradient_per_m = .1;
  p.attitude_balance_rear_gradient_per_m = -.1;
  p.attitude_lift_delta_min = p.attitude_balance_delta_min = -.2;
  p.attitude_lift_delta_max = p.attitude_balance_delta_max = .2;
  const auto datum = model.evaluate(input);
  if (std::abs(datum.force_cog_N.z - plain.force_cog_N.z) > 1e-9 ||
      std::abs(datum.torque_Nm.y - plain.torque_Nm.y) > 1e-9)
    throw std::runtime_error("Attitude map changes its datum");
  input.pitch_angle_rad = .03;
  const auto nose_down = model.evaluate(input);
  input.pitch_angle_rad = -.03;
  const auto nose_up = model.evaluate(input);
  if (std::abs(nose_down.force_cog_N.z - nose_up.force_cog_N.z) < 1. ||
      std::abs(nose_down.torque_Nm.y - nose_up.torque_Nm.y) < 1.)
    throw std::runtime_error("Attitude coupling disconnected");
  for (double speed : {0., 5., 30., 70.})
    for (double pitch : {-.3, 0., .3}) {
      input.vx_mps = speed;
      input.pitch_angle_rad = pitch;
      input.pitch_rate_radps = .2;
      input.heave_velocity_mps = -.1;
      const auto result = model.evaluate(input);
      if (!std::isfinite(result.force_cog_N.z) || !std::isfinite(result.torque_Nm.y))
        throw std::runtime_error("Attitude aero produced nonfinite response");
      if (speed > 0. &&
          !(result.force_cog_N.z >= .8 * .5 * p.c_l * p.air_density_kgpm3 * p.A_m2 * speed * speed &&
            result.force_cog_N.z <= 1.2 * .5 * p.c_l * p.air_density_kgpm3 * p.A_m2 * speed * speed))
        throw std::runtime_error("Attitude aero left bounded force range");
    }
}
