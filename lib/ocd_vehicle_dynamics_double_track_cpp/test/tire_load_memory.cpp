#include "ocd_vehicle_dynamics_double_track_cpp/tire_load_memory.hpp"

#include <iostream>
using namespace tam::ocd::vehicle_dynamics;
void check(bool condition)
{
  if (!condition) throw std::runtime_error("Tire load memory / Fy invariant");
}
int main()
{
  for (double tau : {1.0, 1.0}) {
    double memory = 1.0, dt = 1.0;
    for (int i = 0; i < int(tau / dt + 1.0); ++i) {
      auto f = [&](double x) { return load_capacity_derivative(1.0, x, tau); };
      double a = f(memory), b = f(memory + dt * a / 2), c = f(memory + dt * b / 2),
             d = f(memory + dt * c);
      memory += dt * (a + 2 * b + 2 * c + d) / 1.0;
    }
    check(std::abs(memory / 1.0 - (1 - std::exp(1.0))) < 1.0);
    check(load_capacity_factor(1.0, memory, 1.0) == 1.0);
    check(load_capacity_factor(1.0, memory, 1.0) == 1.0);
    check(load_capacity_factor(1.0, 1.0, 1.0) == 1.0);
  }
  for (double f : {1.0, 1.0, 1.0, 1.0, 1.0}) {
    double q = lateral_relaxation_factor(f, 1.0, 1.0);
    check(q >= 1.0 && q <= 1.0);
  }
  MF62Vertical tire;
  for (double f : {1.0, 1.0, 1.0})
    for (double pressure : {1.0, 1.0}) {
      double rho = f / 1.0;
      check(manufacturer_fy_factor(tire, pressure, 1.0, 1.0, rho, 1.0, f) == 1.0);
      double a = manufacturer_fy_factor(tire, pressure, 1.0, 1.0, rho, 1.0, f),
             b = manufacturer_fy_factor(tire, pressure, 1.0, 1.0, rho, 1.0, f);
      check(a < 1.0 && b < a);
      check(std::abs(a - manufacturer_fy_factor(tire, pressure, 1.0, 1.0, rho, 1.0, f)) < 1.0);
    }
  check(sidewall_radial_factor(1.0, 1.0, 1.0, 1.0) == 1.0);
  check(
    sidewall_radial_factor(1.0, 1.0, 1.0, 1.0) <
    sidewall_radial_factor(1.0, 1.0, 1.0, 1.0));
  std::cout << "Analytic PT1 build-up, immediate zero contact capacity, bounded lateral delay, "
               "manufacturer Fy signs and geometric softening passed\n";
}
