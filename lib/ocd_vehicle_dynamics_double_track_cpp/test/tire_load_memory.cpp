#include <cmath>
#include <stdexcept>

#include "ocd_vehicle_dynamics_double_track_cpp/tire_load_memory.hpp"
using namespace tam::ocd::vehicle_dynamics;

int main()
{
  const double normal = 3000., tau = .02, dt = .0001;
  double memory = 0.;
  for (int i = 0; i < 400; ++i) {
    const auto f = [&](double x) { return load_capacity_derivative(normal, x, tau, tau); };
    const double a = f(memory), b = f(memory + dt * a / 2.),
                 c = f(memory + dt * b / 2.), d = f(memory + dt * c);
    memory += dt * (a + 2. * b + 2. * c + d) / 6.;
  }
  if (std::abs(memory - normal * (1. - std::exp(-.04 / tau))) > 1e-6)
    throw std::runtime_error("Load-memory PT1 buildup mismatch");
  if (load_capacity_factor(0., memory, .5) != 0. ||
      !(load_capacity_factor(normal, memory, .5) > .5 &&
        load_capacity_factor(normal, memory, .5) < 1.))
    throw std::runtime_error("Airborne or recovering contact capacity invalid");
  if (lateral_relaxation_factor(normal, 0., 0.) != 1.)
    throw std::runtime_error("Disabled lateral load sensitivity requires a reference load");
  for (double load : {0., 1000., 3000., 6000.}) {
    const double f = lateral_relaxation_factor(load, normal, .1);
    if (!(f >= .9 && f <= 1.1))
      throw std::runtime_error("Lateral delay leaves bounded range");
  }
  if (sidewall_radial_factor(0., .1, normal, 100000.) != 1. ||
      !(sidewall_radial_factor(.005, .1, normal, 100000.) < 1.))
    throw std::runtime_error("Sidewall geometric feedback sign invalid");
}
