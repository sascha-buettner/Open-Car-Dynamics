#include "ocd_vehicle_dynamics_double_track_cpp/generic_radial_lookup.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/generic_radial_compliance.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

using tam::ocd::vehicle_dynamics::GenericRadialCompliance;
using tam::ocd::vehicle_dynamics::GenericRadialLookup;

void require(bool ok)
{
  if (!ok) throw std::runtime_error("Synthetic radial lookup regression");
}

int main()
{
  GenericRadialLookup lookup;
  require(!lookup.configured());
  bool rejected = false;
  try { (void)lookup.stiffness(false, 1., 1., 1.); }
  catch (const std::logic_error &) { rejected = true; }
  require(rejected);

  // Synthetic monotone grid: these values are unrelated to any vehicle or tire data.
  const std::vector<double> speed{1., 2., 3.}, pressure{1., 2., 3., 4.};
  const std::vector<double> load{1., 2., 3., 4., 5., 6., 7., 8.};
  std::vector<double> front(96), rear(96);
  for (int v = 0; v < 3; ++v)
    for (int p = 0; p < 4; ++p)
      for (int f = 0; f < 8; ++f) {
        const int i = (v * 4 + p) * 8 + f;
        front[i] = 10000. + 100. * v + 10. * p + f;
        rear[i] = 20000. + 200. * v + 20. * p + 2. * f;
      }
  lookup.configure(speed, pressure, load, front, rear);
  require(lookup.configured());
  for (bool is_rear : {false, true})
    for (int v = 0; v < 3; ++v)
      for (int p = 0; p < 4; ++p)
        for (int f = 0; f < 8; ++f) {
          const int i = (v * 4 + p) * 8 + f;
          const double expected = is_rear ? rear[i] : front[i];
          require(std::abs(lookup.stiffness(is_rear, speed[v], pressure[p], load[f]) - expected) < 1e-9);
          require(lookup.stiffness(is_rear, -speed[v], pressure[p], load[f]) == expected);
        }
  const double k = lookup.stiffness(false, 2., 2., 4.);
  require(GenericRadialLookup::soften(k, 0., 1., 1.) == k);
  require(GenericRadialLookup::soften(k, .1, 1., 1.) < k);
  rejected = false;
  try { (void)GenericRadialLookup::soften(k, 1., 1., 1.); }
  catch (const std::invalid_argument &) { rejected = true; }
  require(rejected);

  GenericRadialCompliance compliance;
  compliance.reset(lookup, false, 2., 2.);
  constexpr double force = 3500., step = .01;
  const auto lower = compliance.values(force - step), upper = compliance.values(force + step);
  const double tangent = (upper[0] - lower[0]) / (2. * step);
  require(std::abs(tangent * lookup.stiffness(false, 2., 2., force) - 1.) < 1e-5);
  require(upper[0] > lower[0] && upper[1] > lower[1]);
  require(std::abs((upper[1] - lower[1]) / (2. * step) * lookup.stiffness(false, 2., 2., force) / force - 1.) < 1e-5);
}
