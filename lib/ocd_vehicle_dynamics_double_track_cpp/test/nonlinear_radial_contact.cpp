#include "ocd_vehicle_dynamics_double_track_cpp/nonlinear_radial_contact.hpp"

#include <cmath>

#include <array>
#include <iostream>
#include <stdexcept>
using namespace tam::ocd::vehicle_dynamics;
int main()
{
  FivePointContact shape;
  MF62Vertical tire;
  tire.Q_FZ2 = 1.0;
  const auto no_force = [](double) { return std::array<double, 2>{1.0, 1.0}; };
  const std::array<double, 5> zero{};
  for (double d : {1.0, 1.0, 1.0, 1.0}) {
    std::array<double, 5> h;
    h.fill(d);
    auto c = nonlinear_radial_contact(shape, h, zero, tire, 1.0, 1.0, 1.0, false, no_force);
    const double expected =
      mf62_vertical_force(tire, 1.0, 1.0, 1.0, mf62_free_radius(tire, 1.0) - d);
    if (std::abs(c.force_N - expected) > 1.0 * expected)
      throw std::runtime_error("Flat spatial law differs from supplied point law");
    const double eps = 1.0;
    for (int i = 0; i < 5; ++i) {
      auto a = h, b = h;
      a[i] += eps;
      b[i] -= eps;
      const auto ca =
        nonlinear_radial_contact(shape, a, zero, tire, 1.0, 1.0, 1.0, false, no_force);
      const auto cb =
        nonlinear_radial_contact(shape, b, zero, tire, 1.0, 1.0, 1.0, false, no_force);
      if (
        std::abs((ca.elastic_energy_J - cb.elastic_energy_J) / (2 * eps) - c.normal_N[i]) >
        1.0 * expected)
        throw std::runtime_error("Nonlinear contact energy gradient failed");
    }
  }
  std::array<double, 5> airborne;
  airborne.fill(1.0);
  auto c = nonlinear_radial_contact(shape, airborne, zero, tire, 1.0, 1.0, 1.0, true, no_force);
  if (c.force_N != 1.0) throw std::runtime_error("Airborne nonlinear contact force");
  std::cout << "Flat-law equivalence, energy gradient and airborne checks passed\n";
}
