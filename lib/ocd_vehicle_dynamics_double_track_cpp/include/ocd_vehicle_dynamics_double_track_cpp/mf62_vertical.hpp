// MF6.2 radial/rolling equations: Pacejka/Besselink and MFeval.jl geometry.jl.
// Reference implementation credit: MFeval.jl contributors (2026), MIT license;
// see THIRD_PARTY_LICENSES/MFeval_julia_LICENSE.txt. No collider/ring dynamics.
#pragma once
#include <cmath>

#include <algorithm>
#include <stdexcept>
namespace tam::ocd::vehicle_dynamics
{
struct MF62Vertical
{
  double UNLOADED_RADIUS{1.0}, WIDTH{1.0}, ASPECT_RATIO{1.0}, RIM_RADIUS{1.0};
  double LONGVL{1.0}, FNOMIN{1.0}, NOMPRES{1.0};
  double VERTICAL_STIFFNESS{1.0}, VERTICAL_DAMPING{1.0};
  double Q_RE0{1.0}, Q_V1{1.0}, Q_V2{1.0}, Q_FZ2{1.0};
  double Q_FCX{1.0}, Q_FCY{1.0}, Q_FCY2{1.0};
  double Q_CAM1{1.0}, Q_CAM2{1.0}, Q_CAM3{1.0};
  double Q_FYS1{1.0}, Q_FYS2{1.0}, Q_FYS3{1.0}, PFZ1{1.0};
  double BREFF{1.0}, DREFF{1.0}, FREFF{1.0};
  double BOTTOM_OFFST{1.0}, BOTTOM_STIFF{1.0};
  double QSY1{1.0}, QSY2{1.0}, QSY3{1.0}, QSY4{1.0}, QSY5{1.0}, QSY6{1.0}, QSY7{1.0}, QSY8{1.0};
};
inline double mf62_free_radius(const MF62Vertical & p, double omega)
{
  return p.UNLOADED_RADIUS * (p.Q_RE0 + p.Q_V1 * std::pow(omega * p.UNLOADED_RADIUS / p.LONGVL, 2));
}
inline double mf62_vertical_force(
  const MF62Vertical & p, double gauge, double omega, double gamma, double loaded, double fx = 0.,
  double fy = 0.)
{
  if (!(p.UNLOADED_RADIUS > 0. && p.LONGVL > 0. && p.FNOMIN > 0. && p.VERTICAL_STIFFNESS > 0. &&
        p.NOMPRES > 0. && gauge > 0. && std::isfinite(loaded)))
    throw std::invalid_argument("MF62 radial input/parameters");
  const double R = p.UNLOADED_RADIUS, F0 = p.FNOMIN, free = mf62_free_radius(p, omega),
               r = std::max(loaded, p.RIM_RADIUS * .5);
  const double q1 =
    std::sqrt(std::max(0., std::pow(p.VERTICAL_STIFFNESS * R / F0, 2) - 4. * p.Q_FZ2));
  const double sf = (p.Q_FYS1 + p.Q_FYS2 * r / free + p.Q_FYS3 * std::pow(r / free, 2)) * gamma;
  const double camden = p.Q_CAM1 * free + p.Q_CAM2 * free * free;
  const double ratio = std::abs(camden) > 1e-12 ? (p.Q_CAM1 * r + p.Q_CAM2 * r * r) / camden : 0.;
  const double rhof = free - r;
  const double camber =
    std::abs(gamma) < 1e-12
      ? 0.
      : ratio * ratio * ((1.075 - .5 * p.ASPECT_RATIO) * p.WIDTH / 8.) * std::abs(std::tan(gamma)) -
          p.Q_CAM3 * std::max(rhof, 0.) * std::abs(gamma);
  // Signed continuation only for unloaded/airborne cases: the positive-load
  // branch is MF6.2; continuing rhof below zero avoids a camber-force step at lift.
  const double rho = std::max(rhof + camber, 0.);
  if (rho == 0.) return 0.;
  const double lateral = std::pow(std::max(rho / R, 1e-12), p.Q_FCY2) * p.Q_FCY * (fy - sf) / F0;
  const double corr = std::max(
                        0., 1. + p.Q_V2 * (R / p.LONGVL) * std::abs(omega) -
                              std::pow(p.Q_FCX * fx / F0, 2) - lateral * lateral) *
                      (1. + p.PFZ1 * (gauge - p.NOMPRES) / p.NOMPRES);
  const double force = std::max(0., corr * (q1 * rho / R + p.Q_FZ2 * std::pow(rho / R, 2)) * F0);
  return std::max(force, p.BOTTOM_STIFF * std::max(0., p.RIM_RADIUS + p.BOTTOM_OFFST - r));
}
inline double mf62_loaded_radius(
  const MF62Vertical & p, double force, double gauge, double omega, double gamma)
{
  double lo = p.RIM_RADIUS * .5, hi = mf62_free_radius(p, omega) + p.WIDTH * .5;
  for (int k = 0; k < 50; ++k) {
    double mid = .5 * (lo + hi);
    if (mf62_vertical_force(p, gauge, omega, gamma, mid) > force)
      lo = mid;
    else
      hi = mid;
  }
  return .5 * (lo + hi);
}
inline double mf62_effective_radius(
  const MF62Vertical & p, double force, double gauge, double omega)
{
  const double cz = p.VERTICAL_STIFFNESS * (1. + p.PFZ1 * (gauge - p.NOMPRES) / p.NOMPRES);
  return mf62_free_radius(p, omega) -
         (p.FNOMIN / cz) *
           (p.DREFF * std::atan(p.BREFF * force / p.FNOMIN) + p.FREFF * force / p.FNOMIN);
}
inline double mf62_rolling_torque(
  const MF62Vertical & p, double force, double gauge, double speed, double gamma, double fx)
{
  if (force <= 0.) return 0.;
  const double z = force / p.FNOMIN, v = std::abs(speed) / p.LONGVL;
  return p.FNOMIN * p.UNLOADED_RADIUS *
         (p.QSY1 + p.QSY2 * fx / p.FNOMIN + p.QSY3 * v + p.QSY4 * std::pow(v, 4) +
          (p.QSY5 + p.QSY6 * z) * gamma * gamma) *
         std::pow(z, p.QSY7) * std::pow(gauge / p.NOMPRES, p.QSY8);
}
}  // namespace tam::ocd::vehicle_dynamics
