// Copyright 2026 TUM Autonomous Motorsport
#pragma once
#include <cmath>

#include <array>
#include <stdexcept>
namespace tam::ocd::vehicle_dynamics::mf_schwift
{
// Independent small-displacement rigid-ring kernel in nonrotating hub axes.
// q = (longitudinal, lateral, radial, camber, wind-up, yaw).
// M qdd + (D + G) qd + K q = contact_wrench - M hub_acceleration.
// G is skew symmetric: spin gyroscopic coupling supplies no power.
// Sources: op het Veld 2008, equations 2.13--2.18; Besselink et al. 2010,
// DOI 10.1080/00423111003748088. This is not the vendor MF-Swift kernel.
// The caller owns mass/inertia partition, contact kinematics and reaction paths.
using Six = std::array<double, 6>;
struct Parameters
{
  double mass_kg{1.0};
  double transverse_inertia_kgm2{1.0};
  double spin_inertia_kgm2{1.0};
  double longitudinal_vertical_hz{1.0};
  double lateral_hz{1.0};
  double yaw_camber_hz{1.0};
  double windup_hz{1.0};
  Six damping_ratio{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
};
struct Result
{
  Six acceleration{}, hub_reaction{};
  double stored_energy_J{}, dissipation_W{};
};
inline Result evaluate(
  const Parameters & p, const Six & q, const Six & velocity, const Six & hub_acceleration,
  const Six & contact_wrench, double spin_radps)
{
  const Six mass{p.mass_kg,           p.mass_kg,
                 p.mass_kg,           p.transverse_inertia_kgm2,
                 p.spin_inertia_kgm2, p.transverse_inertia_kgm2};
  const Six frequency{
    p.longitudinal_vertical_hz,
    p.lateral_hz,
    p.longitudinal_vertical_hz,
    p.yaw_camber_hz,
    p.windup_hz,
    p.yaw_camber_hz};
  if (!std::isfinite(spin_radps)) throw std::invalid_argument("Nonfinite ring spin");
  Result out;
  Six gyro{};
  gyro[3] = p.spin_inertia_kgm2 * spin_radps * velocity[5];
  gyro[5] = -p.spin_inertia_kgm2 * spin_radps * velocity[3];
  for (int i = 0; i < 6; ++i) {
    if (!(std::isfinite(mass[i]) && mass[i] > 0. && std::isfinite(frequency[i]) &&
          frequency[i] > 0. && std::isfinite(p.damping_ratio[i]) && p.damping_ratio[i] >= 0.))
      throw std::invalid_argument("Invalid rigid-ring structural parameter");
    for (double value : {q[i], velocity[i], hub_acceleration[i], contact_wrench[i]})
      if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite rigid-ring boundary");
    const double omega = 2. * std::acos(-1.) * frequency[i];
    const double k = mass[i] * omega * omega;
    const double d = 2. * p.damping_ratio[i] * mass[i] * omega;
    const double support = k * q[i] + d * velocity[i];
    out.acceleration[i] = (contact_wrench[i] - support - gyro[i]) / mass[i] - hub_acceleration[i];
    // Equal/opposite sidewall support; G belongs to ring angular momentum.
    out.hub_reaction[i] = support;
    out.stored_energy_J += .5 * (mass[i] * velocity[i] * velocity[i] + k * q[i] * q[i]);
    out.dissipation_W += d * velocity[i] * velocity[i];
  }
  return out;
}
inline double residual_contact_stiffness(double measured_total, double carcass)
{
  if (!(std::isfinite(measured_total) && std::isfinite(carcass) && measured_total > 0. &&
        carcass > measured_total))
    throw std::domain_error("Ring/contact series stiffness requires carcass > measured total > 0");
  return measured_total * carcass / (carcass - measured_total);
}
}  // namespace tam::ocd::vehicle_dynamics::mf_schwift
