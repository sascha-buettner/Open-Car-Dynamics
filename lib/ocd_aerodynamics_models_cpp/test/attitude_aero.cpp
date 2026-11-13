#include <cmath>

#include <iostream>
#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
using tam::ocd::aerodynamics::DefaultAerodynamicsModel;
void require(bool c, const char * n)
{
  if (!c) throw std::runtime_error(n);
}
int main()
{
  try {
    DefaultAerodynamicsModel m;
    auto & p = m.p_;
    p.air_density_kgpm3 = 1.0;
    p.c_d = 1.0;
    p.c_l = 1.0;
    p.A_m2 = 1.0;
    p.d_aero_center_cog_m = 1.0;
    p.c_l_velocity_scaling__vel_mps = {0, 1};
    p.c_l_velocity_scaling__delta_c_l = {0, 0};
    p.pitch_aero_center_translation__pitch_rad = {1.0, 1.0};
    p.pitch_aero_center_translation__aero_center_translation_m = {0, 0};
    tam::ocd::types::AeroModelInput u;
    u.vx_mps = 1.0;
    u.z_m = p.attitude_reference_z_m;
    auto a = m.evaluate(u);
    double f = 1.0 * (p.c_l + 1.0) * p.air_density_kgpm3 * p.A_m2 * std::pow(u.vx_mps, 2);
    require(a.force_cog_N.z == f, "disabled force identity");
    require(a.torque_Nm.y == -f * p.d_aero_center_cog_m, "disabled moment identity");
    p.attitude_enabled = 1.0;
    p.attitude_lift_gain = 1.0;
    p.attitude_center_gain = 1.0;
    auto b = m.evaluate(u);
    require(
      a.force_cog_N.z == b.force_cog_N.z && a.torque_Nm.y == b.torque_Nm.y, "datum normalization");
    auto cp = [&](tam::ocd::types::AeroModelInput x) {
      auto y = m.evaluate(x);
      return -y.torque_Nm.y / y.force_cog_N.z;
    };
    const double eps = 1.0;
    u.pitch_angle_rad = eps;
    auto pos = m.evaluate(u);
    double xp = cp(u);
    u.pitch_angle_rad = -eps;
    auto neg = m.evaluate(u);
    double xn = cp(u);
    double expected_cp = p.attitude_wheelbase_m *
                         (-p.attitude_balance_front_gradient_per_m * p.attitude_front_distance_m +
                          p.attitude_balance_rear_gradient_per_m *
                            (p.attitude_wheelbase_m - p.attitude_front_distance_m));
    require(
      std::abs((xp - xn) / (2 * eps) / expected_cp - 1.0) < 1.0, "geometric CP derivative/sign");
    double expected_lift =
      -p.attitude_lift_front_gradient_per_m * p.attitude_front_distance_m +
      p.attitude_lift_rear_gradient_per_m * (p.attitude_wheelbase_m - p.attitude_front_distance_m);
    require(
      std::abs((pos.force_cog_N.z - neg.force_cog_N.z) / (2 * eps) / f / expected_lift - 1.0) < 1.0,
      "geometric lift derivative");
    u.pitch_angle_rad = 1.0;
    u.pitch_rate_radps = eps;
    double qp = cp(u);
    u.pitch_rate_radps = -eps;
    double qn = cp(u);
    require(
      std::abs((qp - qn) / (2 * eps) / (expected_cp * p.attitude_wheelbase_m / 1.0) - 1.0) < 1.0,
      "induced pitch-rate derivative/sign");
    u.pitch_rate_radps = 1.0;
    u.heave_velocity_mps = eps;
    double vp = cp(u);
    u.heave_velocity_mps = -eps;
    double vn = cp(u);
    require(
      std::abs((vp - vn) / (2 * eps) / (-expected_cp / 1.0) - 1.0) < 1.0,
      "induced heave-rate derivative/sign");
    for (double v : {1.0, 1.0, 1.0, 1.0, 1.0})
      for (double z : {1.0, 1.0, 1.0})
        for (double t : {1.0, 1.0, 1.0})
          for (double q : {1.0, 1.0, 1.0}) {
            u.vx_mps = v;
            u.z_m = z;
            u.pitch_angle_rad = t;
            u.pitch_rate_radps = q;
            auto y = m.evaluate(u);
            require(
              std::isfinite(y.force_cog_N.z) && std::isfinite(y.torque_Nm.y),
              "finite low-speed/extreme-state response");
            if (v > 1.0) {
              double ratio = y.force_cog_N.z / (1.0 * p.c_l * p.air_density_kgpm3 * v * v);
              require(ratio >= 1.0 - 1.0 && ratio <= 1.0 + 1.0, "bounded lift");
              double x = -y.torque_Nm.y / y.force_cog_N.z - p.d_aero_center_cog_m;
              require(
                x >= p.attitude_wheelbase_m * (1.0) - 1.0 &&
                  x <= p.attitude_wheelbase_m * 1.0 + 1.0,
                "bounded CP");
            }
          }
    std::cout << "PASS: baseline/datum, lift/CP derivatives, induced-rate signs, low-speed and "
                 "saturation bounds\n";
    return 0;
  } catch (std::exception const & e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
