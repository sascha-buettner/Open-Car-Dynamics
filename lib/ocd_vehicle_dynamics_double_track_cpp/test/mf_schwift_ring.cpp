#include <cmath>
#include <stdexcept>

#include "ocd_vehicle_dynamics_double_track_cpp/mf_schwift_ring.hpp"
using namespace tam::ocd::vehicle_dynamics::mf_schwift;

int main()
{
  Parameters p;
  p.mass_kg = 5.;
  p.transverse_inertia_kgm2 = .10;
  p.spin_inertia_kgm2 = .05;
  p.longitudinal_vertical_hz = 15.;
  p.lateral_hz = 12.;
  p.yaw_camber_hz = 20.;
  p.windup_hz = 18.;
  p.damping_ratio.fill(.05);
  const Six q{.002, -.003, .001, .01, -.02, .015};
  const Six v{.1, -.2, .05, .2, -.1, .15};
  const Six force{10., -20., 5., 1., -2., 3.}, zero{};
  const double spin = 30.;
  const auto r = evaluate(p, q, v, zero, force, spin);
  const Six mass{p.mass_kg, p.mass_kg, p.mass_kg, p.transverse_inertia_kgm2,
                 p.spin_inertia_kgm2, p.transverse_inertia_kgm2};
  const Six frequency{p.longitudinal_vertical_hz, p.lateral_hz,
                      p.longitudinal_vertical_hz, p.yaw_camber_hz,
                      p.windup_hz, p.yaw_camber_hz};
  double energy_rate = 0.;
  for (int i = 0; i < 6; ++i) {
    const double omega = 2. * std::acos(-1.) * frequency[i];
    energy_rate += mass[i] * v[i] * r.acceleration[i] +
                   mass[i] * omega * omega * q[i] * v[i] - force[i] * v[i];
  }
  if (std::abs(energy_rate + r.dissipation_W) > 1e-9)
    throw std::runtime_error("Ring storage/loss identity failed");
  Six yaw_velocity{};
  yaw_velocity[5] = .2;
  if (evaluate(p, zero, yaw_velocity, zero, zero, spin).acceleration[3] <= 0.)
    throw std::runtime_error("Spin/yaw precession sign incorrect");
  for (double a : evaluate(p, zero, zero, zero, zero, spin).acceleration)
    if (a != 0.) throw std::runtime_error("Spin creates force at rest");
  const double total = 100000., carcass = 200000.;
  const double patch = residual_contact_stiffness(total, carcass);
  if (std::abs(1. / (1. / patch + 1. / carcass) - total) > 1e-8)
    throw std::runtime_error("Series compliance double counted");
  bool rejected = false;
  try {
    residual_contact_stiffness(carcass, carcass);
  } catch (const std::domain_error &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Nonpositive residual stiffness accepted");
}
