// Copyright 2026 TUM Autonomous Motorsport
#pragma once
#include <cmath>

#include <algorithm>
#include <memory>
namespace tam::ocd::types
{
struct RoadSurfaceSample
{
  double height_m{0.0};
  double gradient_x{0.0}, gradient_y{0.0};
  bool valid{false};
  double heading_derivative_m{0.};
  bool bridged_missing_centre{false};
};
// Immutable physical road boundary evaluated synchronously at ODE-stage positions.
// Coordinates and gradients are in the yaw-aligned integration road plane.
// Missing coverage is explicit; callers must not silently extrapolate a cliff.
class RoadSurface
{
public:
  virtual ~RoadSurface() = default;
  virtual void prepare(double, double, double) const {}
  virtual RoadSurfaceSample sample(double x_m, double y_m) const = 0;
  virtual RoadSurfaceSample tire_contact(
    double x, double y, double heading, double radius, double micro_radius) const
  {
    if (radius <= 0.) return footprint(x, y, micro_radius);
    const double c = std::cos(heading), s = std::sin(heading), limit = .995 * radius;
    auto best = sample(x, y);
    if (!best.valid) return {};
    auto evaluate = [&](double u) {
      auto q = sample(x + c * u, y + s * u);
      if (q.valid) {
        q.height_m -= radius - std::sqrt(radius * radius - u * u);
        q.heading_derivative_m = u * (-s * q.gradient_x + c * q.gradient_y);
      }
      return q;
    };
    auto derivative = [&](double u) {
      auto q = sample(x + c * u, y + s * u);
      return q.valid ? c * q.gradient_x + s * q.gradient_y - u / std::sqrt(radius * radius - u * u)
                     : NAN;
    };
    double previous = -limit, df = derivative(previous);
    for (int i = 0; i <= 24; ++i) {
      double u = -limit + 2 * limit * i / 24.;
      auto q = evaluate(u);
      if (q.valid && q.height_m > best.height_m) best = q;
      const double d = derivative(u);
      if (i && std::isfinite(df) && std::isfinite(d) && df > 0. && d < 0.) {
        double lo = previous, hi = u;
        for (int k = 0; k < 18; ++k) {
          double mid = .5 * (lo + hi);
          if (derivative(mid) > 0.)
            lo = mid;
          else
            hi = mid;
        }
        q = evaluate(.5 * (lo + hi));
        if (q.valid && q.height_m > best.height_m) best = q;
      }
      previous = u;
      df = d;
    }
    return best;
  }
  virtual RoadSurfaceSample footprint(double x_m, double y_m, double radius_m) const
  {
    if (radius_m <= 0) return sample(x_m, y_m);
    RoadSurfaceSample out;
    out.valid = true;
    for (int i = -1; i <= 1; ++i)
      for (int j = -1; j <= 1; ++j) {
        const double weight = (i == 0 ? 2. : 1.) * (j == 0 ? 2. : 1.) / 16.;
        const auto q = sample(x_m + i * radius_m, y_m + j * radius_m);
        out.valid = out.valid && q.valid;
        out.height_m += weight * q.height_m;
        out.gradient_x += weight * q.gradient_x;
        out.gradient_y += weight * q.gradient_y;
      }
    return out;
  }
};
}  // namespace tam::ocd::types
