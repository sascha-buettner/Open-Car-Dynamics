#pragma once
#include <array>
#include <limits>

#include "ocd_vehicle_dynamics_double_track_cpp/generic_radial_lookup.hpp"
namespace tam::ocd::vehicle_dynamics
{
// Integrate the configured tangent compliance: rho(F)=integral_0^F 1/K(f) df.
// Replacing K(F) directly in F=K(F)*rho would add an unintended K' feedback.
// Fixed-pressure/speed energy is E(F)=integral_0^F f/K(f) df.
class GenericRadialCompliance
{
  bool rear_ = false;
  double speed_ = std::numeric_limits<double>::quiet_NaN(), pressure_ = 0.;
  std::array<double, 9> compression_{}, energy_{};
  std::array<double, 2> integral(double a, double b) const
  {
    constexpr std::array<double, 8> x{-.9602898564975363, -.7966664774136267, -.5255324099163290,
                                      -.1834346424956498, .1834346424956498,  .5255324099163290,
                                      .7966664774136267,  .9602898564975363};
    constexpr std::array<double, 8> w{.1012285362903763, .2223810344533745, .3137066458778873,
                                      .3626837833783620, .3626837833783620, .3137066458778873,
                                      .2223810344533745, .1012285362903763};
    const double width = (b - a) / 2., half = width / 2.;
    std::array<double, 2> value{};
    for (int panel = 0; panel < 2; ++panel) {
      const double middle = a + (panel + .5) * width;
      for (int i = 0; i < 8; ++i) {
        const double f = middle + half * x[i],
                     c = half * w[i] / GenericRadialLookup::stiffness(rear_, speed_, pressure_, f);
        value[0] += c;
        value[1] += f * c;
      }
    }
    return value;
  }

public:
  void reset(bool rear, double speed, double pressure)
  {
    if (rear == rear_ && speed == speed_ && pressure == pressure_) return;
    rear_ = rear;
    speed_ = speed;
    pressure_ = pressure;
    compression_.fill(0.);
    energy_.fill(0.);
    for (int i = 0; i < 8; ++i) {
      const auto v = integral(1000. * i, 1000. * (i + 1));
      compression_[i + 1] = compression_[i] + v[0];
      energy_[i + 1] = energy_[i] + v[1];
    }
  }
  std::array<double, 2> values(double force) const
  {
    if (!(std::isfinite(force) && force >= 0.))
      throw std::invalid_argument("Radial compliance requires nonnegative normal load");
    const int j = std::min(8, int(force / 1000.));
    const auto tail = integral(1000. * j, force);
    return {compression_[j] + tail[0], energy_[j] + tail[1]};
  }
};
}  // namespace tam::ocd::vehicle_dynamics
