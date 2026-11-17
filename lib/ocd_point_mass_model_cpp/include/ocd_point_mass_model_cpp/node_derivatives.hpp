// Copyright 2026 Sascha Buettner
#pragma once
#include <array>

#include "ocd_point_mass_model_cpp/dual2.hpp"
#include "ocd_point_mass_model_cpp/point_mass_diamond.hpp"

namespace tam::ocd::point_mass
{
/// All outputs of the model at one collocation node with their exact derivatives in
/// z = (x, u): outputs dx/ds (5), q_t, q_reg, path rows (n_path).
struct NodeDerivatives
{
  static constexpr int nz = PointMassDiamond::nz;
  static constexpr int nh = nz * (nz + 1) / 2;
  static constexpr int kMaxOut = PointMassDiamond::nx + 2 + 7;
  static constexpr int kQt = PointMassDiamond::nx;
  static constexpr int kQreg = PointMassDiamond::nx + 1;
  static constexpr int kPath = PointMassDiamond::nx + 2;
  int n_out{0};
  std::array<double, kMaxOut> val{};
  std::array<std::array<double, nz>, kMaxOut> jac{};
  std::array<std::array<double, nh>, kMaxOut> hess{};  // packed upper triangle
};

/// Values only (outputs in the order of NodeDerivatives).
inline int evaluate_node(
  const PointMassDiamond & m, const double * x, const double * u, const TrackPoint & t,
  double * out)
{
  double q_t, q_reg;
  m.dynamics(x, u, t, out, q_t, q_reg);
  out[NodeDerivatives::kQt] = q_t;
  out[NodeDerivatives::kQreg] = q_reg;
  m.path(x, u, t, out + NodeDerivatives::kPath);
  return NodeDerivatives::kPath + m.n_path();
}

/// Values, Jacobian and Hessian of every output in one second-order forward pass.
inline void evaluate_node_derivatives(
  const PointMassDiamond & m, const double * x, const double * u, const TrackPoint & t,
  NodeDerivatives & d)
{
  using D = Dual2<PointMassDiamond::nz>;
  D z[PointMassDiamond::nz];
  for (int i = 0; i < PointMassDiamond::nx; ++i) z[i] = D::variable(x[i], i);
  for (int i = 0; i < PointMassDiamond::nu; ++i)
    z[PointMassDiamond::nx + i] = D::variable(u[i], PointMassDiamond::nx + i);
  D out[NodeDerivatives::kMaxOut];
  m.dynamics(z, z + PointMassDiamond::nx, t, out, out[NodeDerivatives::kQt],
             out[NodeDerivatives::kQreg]);
  m.path(z, z + PointMassDiamond::nx, t, out + NodeDerivatives::kPath);
  d.n_out = NodeDerivatives::kPath + m.n_path();
  for (int k = 0; k < d.n_out; ++k) {
    d.val[k] = out[k].v;
    d.jac[k] = out[k].g;
    d.hess[k] = out[k].h;
  }
}

/// sum_k w_k * Hessian_k (packed), e.g. the node's part of the Hessian of the Lagrangian.
inline void weighted_hessian(
  const NodeDerivatives & d, const double * w, std::array<double, NodeDerivatives::nh> & out)
{
  out.fill(0.0);
  for (int k = 0; k < d.n_out; ++k) {
    if (w[k] == 0.0) continue;
    for (int i = 0; i < NodeDerivatives::nh; ++i) out[i] += w[k] * d.hess[k][i];
  }
}

}  // namespace tam::ocd::point_mass
