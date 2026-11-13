// Copyright 2026 TUM Autonomous Motorsport
#pragma once
#include <cmath>

#include <array>
#include <eigen3/Eigen/Dense>
#include <stdexcept>
namespace tam::ocd::vehicle_dynamics
{
// Five-node, massless condensed elastic belt. Not a rubber material FEM.
// Energy: 1/2 u^T K u, K=diag(area weights)+coupling*graph Laplacian.
// Unilateral road constraints u_i >= road_i - hub + static_deflection.
// Positive definiteness gives a unique, passive elastic contact solution. All
// free-node inverse matrices are precomputed; there is no allocation in evaluate.
class FivePointContact
{
  using Matrix = Eigen::Matrix<double, 5, 5>;
  using Vector = Eigen::Matrix<double, 5, 1>;
  Matrix stiffness_;
  std::array<Matrix, 32> response_{};

public:
  struct Result
  {
    std::array<double, 5> normal_N{};
    std::array<double, 5> belt_displacement_m{};
    double force_N{0.0}, elastic_energy_J{0.0}, dissipated_power_W{0.0};
    unsigned active_mask{0};
  };
  explicit FivePointContact(double coupling = 0.04)
  {
    if (!std::isfinite(coupling) || coupling < 0)
      throw std::invalid_argument("belt coupling must be nonnegative");
    stiffness_ = Matrix::Zero();
    stiffness_.diagonal() << 0.4, 0.15, 0.15, 0.15, 0.15;
    const int edges[8][2] = {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {2, 4}, {4, 3}, {3, 1}};
    for (const auto & e : edges) {
      stiffness_(e[0], e[0]) += coupling;
      stiffness_(e[1], e[1]) += coupling;
      stiffness_(e[0], e[1]) -= coupling;
      stiffness_(e[1], e[0]) -= coupling;
    }
    for (unsigned mask = 0; mask < 32; ++mask) {
      Matrix system = stiffness_, rhs = Matrix::Zero();
      for (int i = 0; i < 5; ++i)
        if (mask & (1u << i)) {
          system.row(i).setZero();
          system(i, i) = 1.;
          rhs(i, i) = 1.;
        }
      response_[mask] = system.partialPivLu().solve(rhs);
    }
  }
  Result evaluate(
    const std::array<double, 5> & compression_m, const std::array<double, 5> & compression_rate_mps,
    double stiffness_Npm, double damping_Nspm) const
  {
    Vector h, rate;
    for (int i = 0; i < 5; ++i) {
      h[i] = compression_m[i];
      rate[i] = compression_rate_mps[i];
    }
    if (
      !h.allFinite() || !rate.allFinite() || !std::isfinite(stiffness_Npm) || stiffness_Npm <= 0 ||
      !std::isfinite(damping_Nspm) || damping_Nspm < 0)
      throw std::invalid_argument("invalid five-point tire input");
    Result out;
    // The common all-contact case returns immediately. Only a curb edge/lift-off
    // requires the bounded (32 sets) complementarity search.
    for (int mask = 31; mask >= 0; --mask) {
      Vector u = response_[mask] * h, force = stiffness_ * u;
      bool feasible = true;
      for (int i = 0; i < 5; ++i) {
        if (mask & (1 << i)) {
          if (force[i] < -1e-10) feasible = false;
        } else if (u[i] < h[i] - 1e-10)
          feasible = false;
      }
      if (!feasible) continue;
      // Local area-weighted dashpots stay dissipative after unilateral force
      // clipping. Clipping a coupled viscosity matrix need not have that property.
      const double weights[5] = {.4, .15, .15, .15, .15};
      for (int i = 0; i < 5; ++i) {
        out.belt_displacement_m[i] = u[i];
        if (mask & (1 << i))
          out.normal_N[i] =
            std::max(0., stiffness_Npm * force[i] + damping_Nspm * weights[i] * rate[i]);
        out.force_N += out.normal_N[i];
        if (mask & (1 << i))
          out.dissipated_power_W +=
            (out.normal_N[i] - stiffness_Npm * std::max(0., force[i])) * rate[i];
      }
      out.elastic_energy_J = .5 * stiffness_Npm * u.dot(force);
      out.active_mask = mask;
      return out;
    }
    throw std::runtime_error("five-point contact complementarity failed");
  }
};
}  // namespace tam::ocd::vehicle_dynamics
