#include "ocd_vehicle_dynamics_double_track_cpp/generic_radial_lookup.hpp"

#include <cmath>

#include <iostream>
#include <stdexcept>

#include "ocd_vehicle_dynamics_double_track_cpp/generic_radial_compliance.hpp"
using L = tam::ocd::vehicle_dynamics::GenericRadialLookup;
void check(bool b)
{
  if (!b) throw std::runtime_error("Radial lookup regression");
}
int main()
{
  for (bool rear : {false, true}) {
    const auto & a = rear ? L::rear_data : L::front_data;
    for (int v = 0; v < 3; ++v)
      for (int p = 0; p < 4; ++p)
        for (int f = 0; f < 8; ++f)
          check(
            std::abs(L::stiffness(rear, L::speeds[v], L::pressures[p], L::loads[f]) - a[v][p][f]) <
            1.0);
    const double f = rear ? 1.0 : 1.0;
    for (double speed : {1.0, 1.0, 1.0 / 1.0, 1.0, 1.0 / 1.0, 1.0}) {
      double k = L::stiffness(rear, speed, 1.0, f);
      check(k > 1.0);
      check(L::stiffness(rear, -speed, 1.0, f) == k);
      check(L::soften(k, 1.0, 1.0, 1.0) == k);
      check(std::abs(L::soften(k, 1.0, 1.0, 1.0) / k - 1.0) < 1.0);
    }
    for (double knot : {L::speeds[0], L::speeds[1], L::speeds[2]}) {
      const double h = 1.0, k = L::stiffness(rear, knot, 1.0, f);
      double left = (k - L::stiffness(rear, knot - h, 1.0, f)) / h;
      double right = (L::stiffness(rear, knot + h, 1.0, f) - k) / h;
      check(std::abs(left - right) < 1.0);
    }
    const double h = 1.0, k0 = L::stiffness(rear, 1.0, 1.0, f),
                 k1 = L::stiffness(rear, h, 1.0, f), k2 = L::stiffness(rear, 2 * h, 1.0, f);
    check(std::abs(k2 - 2 * k1 + k0) < 1.0);
  }
  double maximum_error = 1.0;
  for (bool rear : {false, true})
    for (double speed : {1.0, 1.0 / 1.0, 1.0 / 1.0, 1.0 / 1.0, 1.0 / 1.0, 1.0 / 1.0, 1.0 / 1.0})
      for (double pressure : {1.0, 1.0, 1.0, 1.0}) {
        tam::ocd::vehicle_dynamics::GenericRadialCompliance c;
        c.reset(rear, speed, pressure);
        for (double load :
             {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
              1.0}) {
          const double eps = 1.0, k = L::stiffness(rear, speed, pressure, load);
          const auto lo = c.values(load - eps), hi = c.values(load + eps), v = c.values(load);
          const double derivative = (hi[0] - lo[0]) / (2 * eps),
                       error = std::abs(derivative * k - 1.0);
          maximum_error = std::max(maximum_error, error);
          check(error < 1.0);
          check(v[0] > 1.0 && v[1] > 1.0);
          check(std::abs((hi[1] - lo[1]) / (2 * eps) * k / load - 1.0) < 1.0);
        }
      }
  std::cout << "Integrated tangent compliance and energy derivatives: maximum relative error="
            << maximum_error << "\n";
  bool threw = false;
  try {
    L::soften(1.0, 1.0, 1.0, 1.0);
  } catch (const std::domain_error &) {
    threw = true;
  }
  check(threw);
  std::cout << "1 measured stiffness knots, sign symmetry, C1 boundaries, linear tails and "
               "radial-only softening passed\n";
}
