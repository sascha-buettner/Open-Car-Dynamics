#pragma once
#include <cmath>

#include <algorithm>
#include <array>
#include <stdexcept>
namespace tam::ocd::vehicle_dynamics
{
// Placeholder radial grid; no measured tire data are included.
// Configured speed knots represent steady rolling speeds.
// This is a radial-only empirical law: no tangential state or force is added.
struct GenericRadialLookup
{
  using LoadRow = std::array<double, 8>;
  using PressureRows = std::array<LoadRow, 4>;
  using SpeedRows = std::array<PressureRows, 3>;
  // Deliberately dummy placeholder grid. Physical tire data must be configured privately.
  static constexpr std::array<double, 3> speeds{1.0, 1.0, 1.0};
  static constexpr std::array<double, 4> pressures{1.0, 1.0, 1.0, 1.0};
  static constexpr std::array<double, 8> loads{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  static inline const SpeedRows front_data = [] {
    SpeedRows data{};
    for (auto & pressure_rows : data)
      for (auto & load_row : pressure_rows)
        load_row.fill(1.0);
    return data;
  }();
  static inline const SpeedRows rear_data = front_data;
  static double endpoint(double a, double b)
  {
    double d = (3. * a - b) / 2.;
    if (d * a <= 0.) return 0.;
    if (a * b < 0. && std::abs(d) > 3. * std::abs(a)) return 3. * a;
    return d;
  }
  static double speed_spline(double speed, const std::array<double, 3> & y)
  {
    const double h = speeds[1] - speeds[0], a = (y[1] - y[0]) / h, b = (y[2] - y[1]) / h;
    const std::array<double, 3> d{
      endpoint(a, b), a * b <= 0. ? 0. : 2. * a * b / (a + b), endpoint(b, a)};
    const double x = std::abs(speed);
    if (x <= speeds[0]) return y[0] + d[0] * (x - speeds[0]);
    if (x >= speeds[2]) return y[2] + d[2] * (x - speeds[2]);
    const int j = x < speeds[1] ? 0 : 1;
    const double t = (x - speeds[j]) / h, t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * y[j] + (t3 - 2 * t2 + t) * h * d[j] +
           (-2 * t3 + 3 * t2) * y[j + 1] + (t3 - t2) * h * d[j + 1];
  }
  template <size_t N>
  static int interval(const std::array<double, N> & x, double value)
  {
    return std::clamp(
      int(std::upper_bound(x.begin(), x.end(), value) - x.begin()) - 1, 0, int(N) - 2);
  }
  static double stiffness(bool rear, double speed, double gauge, double load)
  {
    if (!(std::isfinite(speed) && std::isfinite(gauge) && std::isfinite(load) && gauge > 0. &&
          load >= 0.))
      throw std::invalid_argument("Generic radial lookup input");
    const auto & a = rear ? rear_data : front_data;
    const int p = interval(pressures, gauge), n = interval(loads, load);
    const double wp = (gauge - pressures[p]) / (pressures[p + 1] - pressures[p]);
    const double wn = (load - loads[n]) / (loads[n + 1] - loads[n]);
    std::array<double, 3> values{};
    for (int v = 0; v < 3; ++v) {
      const double low = a[v][p][n] + wn * (a[v][p][n + 1] - a[v][p][n]);
      const double high = a[v][p + 1][n] + wn * (a[v][p + 1][n + 1] - a[v][p + 1][n]);
      values[v] = low + wp * (high - low);
    }
    const double k = speed_spline(speed, values);
    if (!(std::isfinite(k) && k > 1000.))
      throw std::domain_error("Nonpositive/excessively soft extrapolated radial stiffness");
    return k;
  }
  static double soften(double k, double fraction, double fx, double reference)
  {
    if (!(fraction >= 0. && fraction <= .3 && reference > 0. && std::isfinite(fx)))
      throw std::invalid_argument("Radial-only force-conditioning parameters");
    const double factor = 1. - fraction * std::pow(fx / reference, 2);
    if (!(factor > .05)) throw std::domain_error("Radial softening outside positive-force domain");
    return k * factor;
  }
};
}  // namespace tam::ocd::vehicle_dynamics
