// Copyright 2026 TUM Autonomous Motorsport
#pragma once

#include <cmath>

#include <algorithm>
#include <array>
#include <boost/math/tools/toms748_solve.hpp>
#include <stdexcept>

#include "ocd_vehicle_dynamics_double_track_cpp/five_point_contact.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/mf62_vertical.hpp"
namespace tam::ocd::vehicle_dynamics
{
// Spatial extension of the supplied radial force law. The condensed elastic
// shape defines rho=sqrt(2*E_unit). With force-coupling disabled, U(rho)=integral
// F_radial drho and its gradient distributes the exact nonlinear force to nodes.
// This is our contact discretization, not a supplied MF-Swift enveloping kernel.
template <class TangentialForce>
FivePointContact::Result nonlinear_radial_contact(
  const FivePointContact & shape, const std::array<double, 5> & compression,
  const std::array<double, 5> & rate, const MF62Vertical & tire, double pressure, double omega,
  double gamma, bool force_coupling, TangentialForce tangential)
{
  auto out = shape.evaluate(compression, rate, 1., 0.);
  const double rho = std::sqrt(2. * out.elastic_energy_J);
  if (rho < 1e-14 || out.force_N < 1e-14) return out;
  const double shape_ratio = out.force_N / rho;
  const double free_radius = mf62_free_radius(tire, omega);
  auto radial = [&](double deflection, const std::array<double, 2> & f) {
    return mf62_vertical_force(tire, pressure, omega, gamma, free_radius - deflection, f[0], f[1]);
  };
  const std::array<double, 2> zero{0., 0.};
  const double uncoupled = shape_ratio * radial(rho, zero);
  double total = uncoupled;
  std::array<double, 2> forces = zero;
  if (force_coupling && uncoupled > 0.) {
    auto residual = [&](double normal) {
      return normal - shape_ratio * radial(rho, tangential(normal));
    };
    double hi = std::max(1., 2. * uncoupled);
    for (int i = 0; residual(hi) < 0. && i < 12; ++i) hi *= 2.;
    const double lo_residual = residual(0.), hi_residual = residual(hi);
    if (!(lo_residual <= 0. && hi_residual >= 0.))
      throw std::domain_error("MF62 contact has no positive-load bracket");
    boost::uintmax_t iterations = 80;
    const auto bracket = boost::math::tools::toms748_solve(
      residual, 0., hi, lo_residual, hi_residual, boost::math::tools::eps_tolerance<double>(44),
      iterations);
    total = .5 * (bracket.first + bracket.second);
    if (std::abs(residual(total)) > 1e-9 * std::max(1., total))
      throw std::runtime_error("MF62 contact root residual exceeds tolerance");
    forces = tangential(total);
  }
  const double scale = total / out.force_N;
  const double weights[5] = {.4, .15, .15, .15, .15};
  out.force_N = 0.;
  out.dissipated_power_W = 0.;
  for (int i = 0; i < 5; ++i) {
    const double elastic = scale * out.normal_N[i];
    out.normal_N[i] = out.active_mask & (1u << i)
                        ? std::max(0., elastic + tire.VERTICAL_DAMPING * weights[i] * rate[i])
                        : 0.;
    out.force_N += out.normal_N[i];
    out.dissipated_power_W += (out.normal_N[i] - elastic) * rate[i];
  }
  // Eight-point Gauss-Legendre integration of the frozen-condition potential.
  // With live Fx/Fy feedback this diagnostic is conditional stored energy,
  // not a proof of passivity for the complete coupled tire.
  constexpr double abscissa[4] = {
    .1834346424956498, .5255324099163290, .7966664774136267, .9602898564975363};
  constexpr double weight[4] = {
    .3626837833783620, .3137066458778873, .2223810344533745, .1012285362903763};
  out.elastic_energy_J = 0.;
  for (int i = 0; i < 4; ++i) {
    out.elastic_energy_J += .5 * rho * weight[i] *
                            (radial(.5 * rho * (1. - abscissa[i]), forces) +
                             radial(.5 * rho * (1. + abscissa[i]), forces));
  }
  return out;
}
}  // namespace tam::ocd::vehicle_dynamics
