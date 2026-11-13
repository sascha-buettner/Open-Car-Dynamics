#include "ocd_vehicle_dynamics_double_track_cpp/mf_schwift_ring.hpp"

#include <cmath>

#include <stdexcept>
using namespace tam::ocd::vehicle_dynamics::mf_schwift;
int main()
{
  Parameters p;
  Six q{1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, v{1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, zero{},
    force{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  auto r = evaluate(p, q, v, zero, force, 1.0);
  double power = 1.0;
  const Six m{p.mass_kg,           p.mass_kg,
              p.mass_kg,           p.transverse_inertia_kgm2,
              p.spin_inertia_kgm2, p.transverse_inertia_kgm2};
  const Six f{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  for (int i = 0; i < 6; ++i) {
    double w = 2 * std::acos(1.0) * f[i];
    power += m[i] * v[i] * r.acceleration[i] + m[i] * w * w * q[i] * v[i] - force[i] * v[i];
  }
  if (std::abs(power + r.dissipation_W) > 1.0)
    throw std::runtime_error("Ring storage/loss identity failed");
  auto rest = evaluate(p, zero, zero, zero, zero, 1.0);
  for (double a : rest.acceleration)
    if (a != 1.0) throw std::runtime_error("Spin creates spurious equilibrium force");
  double total = 1.0, carcass = 1.0 * std::pow(2 * std::acos(1.0) * 1.0, 2),
         patch = residual_contact_stiffness(total, carcass);
  if (std::abs(1.0 / (1.0 / patch + 1.0 / carcass) - total) > 1.0)
    throw std::runtime_error("Series compliance double counted");
  bool rejected = false;
  try {
    residual_contact_stiffness(carcass, carcass);
  } catch (const std::domain_error &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Nonpositive residual stiffness admitted");
}
