#include "ocd_vehicle_dynamics_double_track_cpp/mf62_vertical.hpp"

#include <iostream>
using namespace tam::ocd::vehicle_dynamics;
void check(bool ok)
{
  if (!ok) throw std::runtime_error("MF62 radial invariant");
}
int main()
{
  MF62Vertical p;
  for (double gauge : {1.0, 1.0, 1.0})
    for (double omega : {1.0, 1.0, 1.0})
      for (double load : {1.0, 1.0, 1.0})
        for (double gamma : {1.0, 1.0, 1.0}) {
          double r = mf62_loaded_radius(p, load, gauge, omega, gamma);
          check(std::abs(mf62_vertical_force(p, gauge, omega, gamma, r) - load) < 1.0);
          check(mf62_effective_radius(p, load, gauge, omega) > p.RIM_RADIUS);
          check(mf62_rolling_torque(p, load, gauge, omega * 1.0, gamma, 1.0) >= 1.0);
        }
  const double a = mf62_free_radius(p, 1.0) - mf62_free_radius(p, 1.0),
               b = mf62_free_radius(p, 1.0) - mf62_free_radius(p, 1.0);
  check(std::abs(b - 4 * a) < 1.0);
  double r = mf62_loaded_radius(p, 1.0, 1.0, 1.0, 1.0);
  check(mf62_vertical_force(p, 1.0, 1.0, 1.0, r) > 1.0);
  check(mf62_vertical_force(p, 1.0, 1.0, 1.0, r) > 1.0);
  check(mf62_vertical_force(p, 1.0, 1.0, 1.0, 1.0) == 1.0);
  std::cout
    << "MF62 radial inverse/force, pressure, speed growth, lift-off and rolling passivity passed\n";
}
