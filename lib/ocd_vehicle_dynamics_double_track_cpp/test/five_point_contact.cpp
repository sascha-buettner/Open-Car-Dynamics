#include <array>
#include <cmath>
#include <random>
#include <stdexcept>

#include "ocd_vehicle_dynamics_double_track_cpp/five_point_contact.hpp"
using tam::ocd::vehicle_dynamics::FivePointContact;

int main()
{
  FivePointContact contact;
  constexpr double stiffness = 200000.;
  std::array<double, 5> h{}, rate{};
  h.fill(.01);
  auto flat = contact.evaluate(h, rate, stiffness, 0.);
  if (std::abs(flat.force_N - 2000.) > 1e-7 ||
      std::abs(flat.elastic_energy_J - 10.) > 1e-8 || flat.active_mask != 31)
    throw std::runtime_error("Flat-road contact invariant failed");
  h.fill(0.);
  if (contact.evaluate(h, rate, stiffness, 0.).force_N != 0.)
    throw std::runtime_error("Separated contact transmits force");
  h = {.02, 0., 0., 0., 0.};
  auto curb = contact.evaluate(h, rate, stiffness, 0.);
  if (curb.active_mask == 31 || curb.force_N <= 0.)
    throw std::runtime_error("Raised curb was not enveloped");
  std::mt19937 generator(42);
  std::uniform_real_distribution<double> height(-.005, .02);
  for (int trial = 0; trial < 100; ++trial) {
    for (double & value : h) value = height(generator);
    const auto result = contact.evaluate(h, rate, stiffness, 0.);
    for (int i = 0; i < 5; ++i) {
      if (result.normal_N[i] < 0. || result.belt_displacement_m[i] < h[i] - 1e-9)
        throw std::runtime_error("Unilateral contact condition failed");
      if (result.normal_N[i] > 1e-6 &&
          std::abs(result.belt_displacement_m[i] - h[i]) > 1e-9)
        throw std::runtime_error("Active node does not touch the road");
    }
    for (int i = 0; i < 5; ++i) {
      auto high = h, low = h;
      high[i] += 1e-7;
      low[i] -= 1e-7;
      const double gradient =
        (contact.evaluate(high, rate, stiffness, 0.).elastic_energy_J -
         contact.evaluate(low, rate, stiffness, 0.).elastic_energy_J) / (2e-7);
      if (std::abs(gradient - result.normal_N[i]) > .01)
        throw std::runtime_error("Contact energy gradient mismatch");
    }
  }
  h.fill(.01);
  rate.fill(.2);
  if (contact.evaluate(h, rate, stiffness, 50.).dissipated_power_W <= 0.)
    throw std::runtime_error("Contact damping creates energy");
}
