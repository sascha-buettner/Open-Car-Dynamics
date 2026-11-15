#pragma once
#include <cmath>

#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>
namespace tam::ocd::vehicle_dynamics
{
// Radial-only empirical law. The measured grid must come from private parameters.
struct GenericRadialLookup
{
  using LoadRow = std::array<double, 8>;
  using PressureRows = std::array<LoadRow, 4>;
  using SpeedRows = std::array<PressureRows, 3>;
  std::array<double, 3> speeds{};
  std::array<double, 4> pressures{};
  std::array<double, 8> loads{};
  SpeedRows front_data{}, rear_data{};
  bool configured_{false};
  bool configured() const { return configured_; }
  void configure(
    const std::vector<double> & speed, const std::vector<double> & pressure,
    const std::vector<double> & load, const std::vector<double> & front,
    const std::vector<double> & rear)
  {
    if (speed.size() != speeds.size() || pressure.size() != pressures.size() ||
        load.size() != loads.size() || front.size() != 3 * 4 * 8 || rear.size() != 3 * 4 * 8)
      throw std::invalid_argument("Radial lookup grid is missing or has the wrong shape");
    const auto ascending = [](const auto & values) {
      if (!std::isfinite(values[0]) || values[0] <= 0.) return false;
      for (std::size_t i = 1; i < values.size(); ++i)
        if (!std::isfinite(values[i]) || values[i] <= values[i - 1]) return false;
      return true;
    };
    if (!ascending(speed) || !ascending(pressure) || !ascending(load))
      throw std::invalid_argument("Radial lookup knots must be positive and increasing");
    for (int v = 0; v < 3; ++v)
      for (int p = 0; p < 4; ++p)
        for (int f = 0; f < 8; ++f) {
          const std::size_t index = (v * 4 + p) * 8 + f;
          if (!(std::isfinite(front[index]) && front[index] > 1000. &&
                std::isfinite(rear[index]) && rear[index] > 1000.))
            throw std::invalid_argument("Radial lookup stiffness must be finite and positive");
          front_data[v][p][f] = front[index];
          rear_data[v][p][f] = rear[index];
        }
    std::copy(speed.begin(), speed.end(), speeds.begin());
    std::copy(pressure.begin(), pressure.end(), pressures.begin());
    std::copy(load.begin(), load.end(), loads.begin());
    configured_ = true;
  }
  static double endpoint(double a, double b)
  {
    double d = (3. * a - b) / 2.;
    if (d * a <= 0.) return 0.;
    if (a * b < 0. && std::abs(d) > 3. * std::abs(a)) return 3. * a;
    return d;
  }
  double speed_spline(double speed, const std::array<double, 3> & y) const
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
  double stiffness(bool rear, double speed, double gauge, double load) const
  {
    if (!configured_) throw std::logic_error("Radial lookup grid was not configured");
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
