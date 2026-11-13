#pragma once
#include <cmath>

#include <algorithm>
#include <stdexcept>

#include "ocd_vehicle_dynamics_double_track_cpp/mf62_vertical.hpp"
namespace tam::ocd::vehicle_dynamics
{
inline double load_capacity_factor(double normal, double remembered, double fraction)
{
  if (!(normal >= 0. && std::isfinite(normal) && std::isfinite(remembered) && fraction >= 0. &&
        fraction <= 1.))
    throw std::domain_error("Invalid load capacity state");
  if (normal <= 0.) return 0.;
  return 1. - fraction + fraction * std::clamp(remembered / normal, 0., 1.);
}
inline double load_capacity_derivative(double normal, double remembered, double tau)
{
  if (!(tau > 0. && std::isfinite(tau)))
    throw std::invalid_argument("Positive load build-up time required");
  return (normal - remembered) / (normal > remembered ? tau : 1.0);
}
inline double lateral_relaxation_factor(double normal, double reference, double strength)
{
  if (!(strength >= 0. && strength <= .1 && reference > 0.))
    throw std::invalid_argument("Small lateral load-memory sensitivity required");
  return 1. + strength * (std::max(0., normal) - reference) / (std::max(0., normal) + reference);
}
// Small two-strut geometric hypothesis, including compression geometric stiffness.
// It is a sensitivity, not a measured carcass constitutive law.
inline double sidewall_radial_factor(double lateral, double height, double normal, double stiffness)
{
  const double r = std::clamp(lateral / height, -.6, .6), c2 = 1. - r * r;
  return c2 - normal / (stiffness * height) * r * r / c2;
}
// Preserve the existing zero-Fy law. Evaluate the manufacturer's Fy/camber
// correction at the same radial deflection, pressure and rolling speed.
inline double manufacturer_fy_factor(
  const MF62Vertical & p, double pressure, double speed, double gamma, double rho, double fy,
  double normal, bool enforce_domain = true)
{
  if (normal < 1e-8) return 1.;
  const double omega = std::abs(speed) / p.UNLOADED_RADIUS;
  const double free = mf62_free_radius(p, omega), loaded = free - std::max(rho, 1e-8);
  const double base = mf62_vertical_force(p, pressure, omega, gamma, loaded, 0., 0.);
  const double changed = mf62_vertical_force(p, pressure, omega, gamma, loaded, 0., fy);
  if (base <= 1e-8) {
    if (normal < 100.) return 1.;
    throw std::domain_error("Undefined zero-Fy manufacturer correction");
  }
  // Fade this correction below the file's declared100N validity boundary.
  const double ratio = 1. + std::min(1., normal / 100.) * (changed / base - 1.);
  if (!std::isfinite(ratio) || (enforce_domain && !(ratio > .05 && ratio < 5.)))
    throw std::domain_error("Manufacturer Fy correction leaves positive radial domain");
  return ratio;
}
}  // namespace tam::ocd::vehicle_dynamics
