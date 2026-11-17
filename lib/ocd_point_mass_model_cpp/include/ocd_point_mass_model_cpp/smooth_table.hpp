// Copyright 2026 Sascha Buettner
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace tam::ocd::point_mass
{
/// Largest number of knots of one table axis (the basis is evaluated into fixed-size buffers).
constexpr int kMaxKnots = 512;

/// Softplus ramp eps * log(1 + exp(t / eps)) and its first two derivatives in t (overflow-safe):
/// it approaches relu(t) exponentially, so the rounding of a knot stays local.
inline void smooth_ramp(double t, double eps, double & r, double & r1, double & r2)
{
  const double u = t / eps;
  const double e = std::exp(-std::fabs(u));
  r = eps * (std::max(u, 0.0) + std::log1p(e));
  r1 = u > 0.0 ? 1.0 / (1.0 + e) : e / (1.0 + e);
  r2 = e / ((1.0 + e) * (1.0 + e)) / eps;
}

/// One axis of a C-infinity version of a piecewise-linear table with linear extrapolation:
/// relu(x - x_k) at each interior knot replaced by the softplus ramp of width
/// eps_k = frac * min(h_{k-1}, h_k). The hat functions phi_i(x) of the table are
///
///     phi_i(x) = e_i(0) + s_i(0) (x - x_0) + sum_k (s_i(k) - s_i(k-1)) ramp_k(x)
///
/// with s_i the segment slopes of the unit table e_i (nonzero at k = i-1, i, i+1 only).
class SmoothAxis
{
public:
  SmoothAxis() = default;
  SmoothAxis(std::vector<double> grid, double frac) : grid_(std::move(grid))
  {
    const int n = size();
    if (n < 1 || n > kMaxKnots) throw std::invalid_argument("SmoothAxis: 1..kMaxKnots knots");
    for (int k = 0; k + 1 < n; ++k)
      if (!(grid_[k + 1] > grid_[k])) throw std::invalid_argument("SmoothAxis: grid must increase");
    eps_.assign(std::max(n - 2, 0), 0.0);
    for (int k = 1; k + 1 < n; ++k)
      eps_[k - 1] = frac * std::min(grid_[k] - grid_[k - 1], grid_[k + 1] - grid_[k]);
  }
  int size() const { return static_cast<int>(grid_.size()); }
  const std::vector<double> & grid() const { return grid_; }

  /// phi_i(x) and, if d1 / d2 are not null, the first and second derivatives (arrays of size()).
  void basis(double x, double * phi, double * d1 = nullptr, double * d2 = nullptr) const
  {
    const int n = size();
    if (n == 1) {
      phi[0] = 1.0;
      if (d1) d1[0] = 0.0;
      if (d2) d2[0] = 0.0;
      return;
    }
    std::array<double, kMaxKnots> r, r1, r2;  // ramps at the interior knots 1 .. n-2
    for (int k = 1; k + 1 < n; ++k) smooth_ramp(x - grid_[k], eps_[k - 1], r[k], r1[k], r2[k]);
    const double dx0 = x - grid_[0];
    for (int i = 0; i < n; ++i) {
      // slopes of the unit table e_i on the segments i-1 and i
      const double s_left = i >= 1 ? 1.0 / (grid_[i] - grid_[i - 1]) : 0.0;
      const double s_right = i + 1 < n ? -1.0 / (grid_[i + 1] - grid_[i]) : 0.0;
      auto slope = [&](int k) { return k == i - 1 ? s_left : (k == i ? s_right : 0.0); };
      const double e0 = i == 0 ? 1.0 : 0.0;
      const double s0 = slope(0);
      double f = e0 + s0 * dx0, f1 = s0, f2 = 0.0;
      for (int k = std::max(1, i - 1); k <= std::min(n - 2, i + 1); ++k) {
        const double ds = slope(k) - slope(k - 1);
        if (ds != 0.0) {
          f += ds * r[k];
          f1 += ds * r1[k];
          f2 += ds * r2[k];
        }
      }
      phi[i] = f;
      if (d1) d1[i] = f1;
      if (d2) d2[i] = f2;
    }
  }

private:
  std::vector<double> grid_;
  std::vector<double> eps_;
};

/// Smooth 1D table v(x) = sum_i v_i phi_i(x).
class SmoothTable1D
{
public:
  SmoothTable1D() = default;
  SmoothTable1D(std::vector<double> grid, std::vector<double> values, double frac)
  : axis_(std::move(grid), frac), values_(std::move(values))
  {
    if (static_cast<int>(values_.size()) != axis_.size())
      throw std::invalid_argument("SmoothTable1D: one value per knot");
  }
  /// Value, first and second derivative at x.
  void eval(double x, double & f, double & f1, double & f2) const
  {
    std::array<double, kMaxKnots> phi, d1, d2;
    axis_.basis(x, phi.data(), d1.data(), d2.data());
    f = f1 = f2 = 0.0;
    for (int i = 0; i < axis_.size(); ++i) {
      f += values_[i] * phi[i];
      f1 += values_[i] * d1[i];
      f2 += values_[i] * d2[i];
    }
  }
  double operator()(double x) const
  {
    std::array<double, kMaxKnots> phi;
    axis_.basis(x, phi.data());
    double f = 0.0;
    for (int i = 0; i < axis_.size(); ++i) f += values_[i] * phi[i];
    return f;
  }
  const SmoothAxis & axis() const { return axis_; }
  const std::vector<double> & values() const { return values_; }

private:
  SmoothAxis axis_;
  std::vector<double> values_;
};

/// Smooth 2D tables on one grid, z(x, y) = phi(x)^T Z psi(y) (tensor product of the two smoothed
/// axes); several tables share the basis evaluation. Z is row-major (x index major).
class SmoothTables2D
{
public:
  struct Derivatives
  {
    double f, fx, fy, fxx, fxy, fyy;
  };

  SmoothTables2D() = default;
  SmoothTables2D(
    std::vector<double> grid_x, std::vector<double> grid_y, std::vector<std::vector<double>> tables,
    double frac)
  : ax_(std::move(grid_x), frac), ay_(std::move(grid_y), frac), tables_(std::move(tables))
  {
    for (const auto & z : tables_)
      if (static_cast<int>(z.size()) != ax_.size() * ay_.size())
        throw std::invalid_argument("SmoothTables2D: table size must be nx * ny");
  }
  int count() const { return static_cast<int>(tables_.size()); }

  /// Values of the tables `which` (indices) at (x, y).
  template <std::size_t K>
  void values(double x, double y, const std::array<int, K> & which, std::array<double, K> & out) const
  {
    std::array<double, kMaxKnots> phi, psi;
    ax_.basis(x, phi.data());
    ay_.basis(y, psi.data());
    const int nx = ax_.size(), ny = ay_.size();
    for (std::size_t t = 0; t < K; ++t) {
      const auto & z = tables_[which[t]];
      double f = 0.0;
      for (int i = 0; i < nx; ++i) {
        double zpsi = 0.0;
        for (int j = 0; j < ny; ++j) zpsi += z[i * ny + j] * psi[j];
        f += phi[i] * zpsi;
      }
      out[t] = f;
    }
  }

  /// Values and derivatives up to second order of the tables `which` at (x, y).
  template <std::size_t K>
  void derivatives(
    double x, double y, const std::array<int, K> & which, std::array<Derivatives, K> & out) const
  {
    std::array<double, kMaxKnots> phi, phi1, phi2, psi, psi1, psi2;
    ax_.basis(x, phi.data(), phi1.data(), phi2.data());
    ay_.basis(y, psi.data(), psi1.data(), psi2.data());
    const int nx = ax_.size(), ny = ay_.size();
    for (std::size_t t = 0; t < K; ++t) {
      const auto & z = tables_[which[t]];
      Derivatives d{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      for (int i = 0; i < nx; ++i) {
        double z0 = 0.0, z1 = 0.0, z2 = 0.0;
        for (int j = 0; j < ny; ++j) {
          const double zij = z[i * ny + j];
          z0 += zij * psi[j];
          z1 += zij * psi1[j];
          z2 += zij * psi2[j];
        }
        d.f += phi[i] * z0;
        d.fx += phi1[i] * z0;
        d.fy += phi[i] * z1;
        d.fxx += phi2[i] * z0;
        d.fxy += phi1[i] * z1;
        d.fyy += phi[i] * z2;
      }
      out[t] = d;
    }
  }

private:
  SmoothAxis ax_, ay_;
  std::vector<std::vector<double>> tables_;
};

}  // namespace tam::ocd::point_mass
