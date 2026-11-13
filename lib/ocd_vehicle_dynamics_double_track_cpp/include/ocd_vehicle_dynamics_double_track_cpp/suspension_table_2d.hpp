// Coupled jounce/rack map: monotone Hermite in jounce, linear in rack.
#pragma once
#include <cmath>

#include <algorithm>
#include <stdexcept>
#include <vector>
namespace tam::ocd::vehicle_dynamics
{
struct SuspensionMap2D
{
  std::vector<double> control{1.0, 1.0};
  std::vector<double> arb{1.0, 1.0, 1.0, 1.0};
  std::vector<double> spring{1.0, 1.0, 1.0, 1.0}, damper{1.0, 1.0, 1.0, 1.0};
  std::vector<double> normal_x{1.0, 1.0, 1.0, 1.0}, normal_y{1.0, 1.0, 1.0, 1.0}, normal_z{1.0, 1.0, 1.0, 1.0};
  std::vector<double> contact_x{1.0, 1.0, 1.0, 1.0}, contact_y{1.0, 1.0, 1.0, 1.0},
    contact_z{1.0, 1.0, 1.0, 1.0};
  std::vector<double> jac_x{1.0, 1.0, 1.0, 1.0}, jac_y{1.0, 1.0, 1.0, 1.0}, jac_z{1.0, 1.0, 1.0, 1.0};
  std::vector<double> angjac_x{1.0, 1.0, 1.0, 1.0}, angjac_y{1.0, 1.0, 1.0, 1.0}, angjac_z{1.0, 1.0, 1.0, 1.0},
    rocker{1.0, 1.0, 1.0, 1.0};
};
struct MapValue
{
  double value, dq, ds;
};
inline std::pair<double, double> map_column(
  double q, const std::vector<double> & x, const std::vector<double> & y, size_t offset)
{
  auto sec = [&](size_t i) { return (y[offset + i + 1] - y[offset + i]) / (x[i + 1] - x[i]); };
  auto tangent = [&](size_t i) {
    if (i == 0) return sec(0);
    if (i + 1 == x.size()) return sec(i - 1);
    const double a = sec(i - 1), b = sec(i);
    return a * b <= 0. ? 0. : 2 * a * b / (a + b);
  };
  if (q < x.front()) return {y[offset] + tangent(0) * (q - x.front()), tangent(0)};
  if (q > x.back())
    return {
      y[offset + x.size() - 1] + tangent(x.size() - 1) * (q - x.back()), tangent(x.size() - 1)};
  const size_t h = std::clamp<size_t>(
                 std::upper_bound(x.begin(), x.end(), q) - x.begin(), 1, x.size() - 1),
               l = h - 1;
  const double dx = x[h] - x[l], t = (q - x[l]) / dx, m0 = tangent(l), m1 = tangent(h),
               yl = y[offset + l], yh = y[offset + h];
  return {
    (2 * t * t * t - 3 * t * t + 1) * yl + (t * t * t - 2 * t * t + t) * dx * m0 +
      (-2 * t * t * t + 3 * t * t) * yh + (t * t * t - t * t) * dx * m1,
    ((6 * t * t - 6 * t) * yl + (3 * t * t - 4 * t + 1) * dx * m0 + (-6 * t * t + 6 * t) * yh +
     (3 * t * t - 2 * t) * dx * m1) /
      dx};
}
inline MapValue suspension_map(
  double q, double s, const std::vector<double> & x, const std::vector<double> & ss,
  const std::vector<double> & y)
{
  if (
    x.size() < 2 || ss.size() < 2 || y.size() != x.size() * ss.size() || !std::isfinite(q) ||
    !std::isfinite(s))
    throw std::invalid_argument("Invalid coupled suspension map");
  const size_t hi = std::clamp<size_t>(
                 std::upper_bound(ss.begin(), ss.end(), s) - ss.begin(), 1, ss.size() - 1),
               lo = hi - 1;
  if (!(ss[hi] > ss[lo])) throw std::invalid_argument("Invalid rack knots");
  const double ds = ss[hi] - ss[lo], t = std::clamp((s - ss[lo]) / ds, 0., 1.);
  const auto a = map_column(q, x, y, lo * x.size()), b = map_column(q, x, y, hi * x.size());
  return {
    a.first + t * (b.first - a.first), a.second + t * (b.second - a.second),
    (s < ss.front() || s > ss.back()) ? 0. : (b.first - a.first) / ds};
}
}  // namespace tam::ocd::vehicle_dynamics
