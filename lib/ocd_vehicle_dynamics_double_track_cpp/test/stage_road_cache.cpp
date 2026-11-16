#include <cmath>
#include <memory>
#include <stdexcept>

#include "ocd_vehicle_dynamics_double_track_cpp/stage_road_cache.hpp"

namespace
{
struct CountingSurface : tam::ocd::types::RoadSurface
{
  mutable int footprint_calls{0};
  mutable int prepare_calls{0};

  void prepare(double, double, double) const override { ++prepare_calls; }

  tam::ocd::types::RoadSurfaceSample sample(double x, double y) const override
  {
    return {x + y, x, y, true};
  }

  tam::ocd::types::RoadSurfaceSample footprint(double x, double y, double radius) const override
  {
    ++footprint_calls;
    return {x + y + radius, x, y, true};
  }
};

void require(bool condition)
{
  if (!condition) throw std::runtime_error("Stage road cache changed exact query behavior");
}
}  // namespace

int main()
{
  auto source = std::make_shared<CountingSurface>();
  tam::ocd::vehicle_dynamics::StageRoadCache cache(source);
  for (int repeat = 0; repeat < 8; ++repeat) {
    cache.prepare(0., 0., 0.);
    for (int i = 0; i < 20; ++i) {
      const auto value = cache.footprint(double(i), 1., .1);
      require(value.valid && std::abs(value.height_m - (i + 1.1)) < 1e-12);
    }
  }
  require(source->prepare_calls == 1 && source->footprint_calls == 20);

  cache.prepare(1., 0., 0.);
  for (int i = 0; i < 20; ++i) cache.footprint(double(i), 1., .1);
  require(source->prepare_calls == 2 && source->footprint_calls == 40);

  cache.prepare(1., 0., 0.);
  for (int i = 0; i < 40; ++i) cache.footprint(double(i), 1., .1);
  require(source->prepare_calls == 2 && source->footprint_calls == 60);
}
