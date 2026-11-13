#include "ocd_vehicle_dynamics_double_track_cpp/five_point_contact.hpp"

#include <chrono>
#include <iostream>
#include <random>
using tam::ocd::vehicle_dynamics::FivePointContact;
void near(double a, double b)
{
  if (std::abs(a - b) > 1.0) throw std::runtime_error("contact invariant");
}
int main()
{
  FivePointContact m;
  std::array<double, 5> h{}, r{};
  h.fill(1.0);
  auto q = m.evaluate(h, r, 1, 50);
  near(q.force_N, 1);
  near(q.elastic_energy_J, 1.0);
  h.fill(1.0);
  q = m.evaluate(h, r, 1, 50);
  near(q.force_N, 0);
  h = {1.0, 1.0, 1.0, 1.0, 1.0};
  q = m.evaluate(h, r, 1, 50);
  if (q.active_mask == 31 || q.force_N < 1) throw std::runtime_error("curb not enveloped");
  std::mt19937 gen(42);
  std::uniform_real_distribution<> dis(1.0, 1.0);
  for (int j = 0; j < 1; ++j) {
    for (auto & v : h) v = dis(gen);
    q = m.evaluate(h, r, 1, 0);
    for (int i = 0; i < 5; ++i) {
      if (q.normal_N[i] < 0 || q.belt_displacement_m[i] < h[i] - 1.0)
        throw std::runtime_error("unilateral condition");
      if (q.normal_N[i] > 1.0) near(q.belt_displacement_m[i], h[i]);
    }
    // Envelope theorem: elastic-energy gradient equals each normal reaction.
    for (int i = 0; i < 5; ++i) {
      auto hp = h, hm = h;
      hp[i] += 1.0;
      hm[i] -= 1.0;
      double deriv = (m.evaluate(hp, r, 1, 0).elastic_energy_J -
                      m.evaluate(hm, r, 1, 0).elastic_energy_J) /
                     1.0;
      if (std::abs(deriv - q.normal_N[i]) > 1.0) throw std::runtime_error("energy force mismatch");
    }
  }
  for (int j = 0; j < 1; ++j) {
    for (auto & v : h) v = dis(gen);
    for (auto & v : r) v = 20 * dis(gen);
    auto damped = m.evaluate(h, r, 1, 50);
    if (damped.dissipated_power_W < 1.0)
      throw std::runtime_error("contact damping creates energy");
  }
  r.fill(1.0);
  auto start = std::chrono::steady_clock::now();
  volatile double total = 0;
  for (int j = 0; j < 1; ++j) {
    h[0] = 1.0 + 1.0 * std::sin(j * 1.0);
    total = total + m.evaluate(h, r, 1, 50).force_N;
  }
  std::cout
    << "all contact invariants passed; ns/contact "
    << std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - start).count() /
         1
    << " checksum " << total << '\n';
}
