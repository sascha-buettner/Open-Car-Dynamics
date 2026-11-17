// Copyright 2026 Sascha Buettner
// Exact derivatives of the point mass against central finite differences of its values.
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "ocd_point_mass_model_cpp/node_derivatives.hpp"

using namespace tam::ocd::point_mass;  // NOLINT

namespace
{
PointMassDiamond make_model(bool yaw)
{
  PointMassParams p;
  p.h_cog_m = 0.3;
  p.Iz_kgm2 = 1000.0;
  p.w_jerk_x = 1.0e-3;
  p.w_jerk_y = 2.0e-3;
  p.w_kappa = 0.1;
  p.gg_exponent_ax_pos = 1.8;
  p.gg_exponent_ax_neg = 1.5;
  p.yaw_moment_limit = yaw;
  std::vector<double> v_list{10, 20, 30, 40, 50, 60, 70, 80};
  std::vector<double> g_list{2, 4, 6, 8, 10, 12, 14};
  std::vector<std::vector<double>> tables(4, std::vector<double>(v_list.size() * g_list.size()));
  for (size_t i = 0; i < v_list.size(); ++i)
    for (size_t j = 0; j < g_list.size(); ++j) {
      const double v = v_list[i], g = g_list[j];
      tables[0][i * g_list.size() + j] = -(1.0 + 0.002 * v * v / 10.0) * g;
      tables[1][i * g_list.size() + j] = 0.9 * g - 0.01 * v;
      tables[2][i * g_list.size() + j] = (1.1 + 0.001 * v) * g;
      tables[3][i * g_list.size() + j] = 3000.0 + 20.0 * g + 5.0 * v;
    }
  DiamondLimits lim{
    SmoothTables2D(v_list, g_list, tables, 0.05),
    SmoothTable1D({0, 26.6, 26.65, 36, 36.05, 90}, {12, 13, 12.2, 9, 8.3, 0}, 0.05)};
  return PointMassDiamond(p, lim);
}
}  // namespace

int main()
{
  int failures = 0;
  std::mt19937 rng(1);
  std::uniform_real_distribution<double> uni(-1.0, 1.0);
  for (bool yaw : {false, true}) {
    const PointMassDiamond m = make_model(yaw);
    for (int trial = 0; trial < 50; ++trial) {
      double x[5] = {35.0 + 20.0 * uni(rng), 3.0 * uni(rng), 0.2 * uni(rng), 8.0 * uni(rng),
                     12.0 * uni(rng)};
      double u[2] = {20.0 * uni(rng), 30.0 * uni(rng)};
      TrackPoint t{0.03 * uni(rng), 0.05 * uni(rng), 1e-3 * uni(rng), 1e-5 * uni(rng),
                   1e-3 * uni(rng), 1e-5 * uni(rng), 0.01 * uni(rng), 1e-4 * uni(rng)};
      NodeDerivatives d;
      evaluate_node_derivatives(m, x, u, t, d);
      double z[7] = {x[0], x[1], x[2], x[3], x[4], u[0], u[1]};
      auto f = [&](const double * zz, double * out) { return evaluate_node(m, zz, zz + 5, t, out); };
      double base[NodeDerivatives::kMaxOut];
      f(z, base);
      double max_val = 0.0, max_jac = 0.0, max_hess = 0.0;
      for (int k = 0; k < d.n_out; ++k)
        max_val = std::max(max_val, std::fabs(base[k] - d.val[k]) / (1.0 + std::fabs(base[k])));
      for (int i = 0; i < 7; ++i) {
        const double hi = 1e-7 * (1.0 + std::fabs(z[i]));  // the 2.5 mm machine step needs a small step
        double zp[7], zm[7], fp[NodeDerivatives::kMaxOut], fm[NodeDerivatives::kMaxOut];
        std::copy(z, z + 7, zp);
        std::copy(z, z + 7, zm);
        zp[i] += hi;
        zm[i] -= hi;
        f(zp, fp);
        f(zm, fm);
        for (int k = 0; k < d.n_out; ++k) {
          const double fd = (fp[k] - fm[k]) / (2.0 * hi);
          max_jac = std::max(max_jac, std::fabs(fd - d.jac[k][i]) / (1.0 + std::fabs(fd)));
        }
        // Hessian column i: central difference of the exact gradient
        NodeDerivatives dp, dm;
        evaluate_node_derivatives(m, zp, zp + 5, t, dp);
        evaluate_node_derivatives(m, zm, zm + 5, t, dm);
        for (int k = 0; k < d.n_out; ++k)
          for (int j = 0; j < 7; ++j) {
            const double fd = (dp.jac[k][j] - dm.jac[k][j]) / (2.0 * hi);
            const int a = std::min(i, j), b = std::max(i, j);
            const double ex = d.hess[k][packed_index(7, a, b)];
            max_hess = std::max(max_hess, std::fabs(fd - ex) / (1.0 + std::fabs(fd)));
          }
      }
      if (max_val > 1e-13 || max_jac > 1e-5 || max_hess > 1e-5) {
        ++failures;
        std::printf("yaw %d trial %d: val %.2e jac %.2e hess %.2e\n", yaw, trial, max_val,
                    max_jac, max_hess);
      }
    }
  }
  std::printf("%s (%d failures)\n", failures ? "FAILED" : "ok", failures);
  return failures ? 1 : 0;
}
