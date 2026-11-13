#pragma once
#include <cmath>

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>
namespace tam::ocd::vehicle_dynamics
{
inline std::pair<double, double> suspension_lookup(
  double q, const std::vector<double> & x, const std::vector<double> & y)
{
  if (x.size() < 2 || x.size() != y.size() || !std::isfinite(q))
    throw std::invalid_argument("Invalid suspension lookup dimensions/input");
  for (size_t i = 0; i < x.size(); ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || (i && x[i] <= x[i - 1]))
      throw std::invalid_argument("Invalid suspension lookup knots");
  if (q < x.front()) return {y.front(), 0.};
  if (q > x.back()) return {y.back(), 0.};
  const size_t h =
    std::clamp<size_t>(std::upper_bound(x.begin(), x.end(), q) - x.begin(), 1, x.size() - 1);
  const double slope = (y[h] - y[h - 1]) / (x[h] - x[h - 1]);
  return {y[h - 1] + slope * (q - x[h - 1]), slope};
}
// Monotone cubic Hermite map, with continuous analytic tangent for virtual-work forces.
inline std::pair<double, double> suspension_smooth_lookup(
  double q, const std::vector<double> & x, const std::vector<double> & y)
{
  suspension_lookup(q, x, y);  // Validate before indexing.
  auto secant = [&](size_t i) { return (y[i + 1] - y[i]) / (x[i + 1] - x[i]); };
  auto tangent = [&](size_t i) {
    if (i == 0) return secant(0);
    if (i + 1 == x.size()) return secant(i - 1);
    const double a = secant(i - 1), b = secant(i);
    return a * b <= 0. ? 0. : 2. * a * b / (a + b);
  };
  if (q < x.front()) return {y.front() + tangent(0) * (q - x.front()), tangent(0)};
  if (q > x.back())
    return {y.back() + tangent(x.size() - 1) * (q - x.back()), tangent(x.size() - 1)};
  const size_t h = std::clamp<size_t>(
                 std::upper_bound(x.begin(), x.end(), q) - x.begin(), 1, x.size() - 1),
               l = h - 1;
  const double dx = x[h] - x[l], t = (q - x[l]) / dx, m0 = tangent(l), m1 = tangent(h);
  return {
    (2 * t * t * t - 3 * t * t + 1) * y[l] + (t * t * t - 2 * t * t + t) * dx * m0 +
      (-2 * t * t * t + 3 * t * t) * y[h] + (t * t * t - t * t) * dx * m1,
    ((6 * t * t - 6 * t) * y[l] + (3 * t * t - 4 * t + 1) * dx * m0 + (-6 * t * t + 6 * t) * y[h] +
     (3 * t * t - 2 * t) * dx * m1) /
      dx};
}
}  // namespace tam::ocd::vehicle_dynamics
