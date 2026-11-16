#pragma once
#include <array>
#include <map>
#include <memory>

#include "ocd_types_cpp/road_surface.hpp"
namespace tam::ocd::vehicle_dynamics
{
// Exact memoization only within one ODE evaluation; no spatial interpolation,
// previous-stage force, or previous-stage surface sample is reused.
class StageRoadCache final : public types::RoadSurface
{
  std::shared_ptr<const types::RoadSurface> source_;
  mutable std::map<std::array<double, 2>, types::RoadSurfaceSample> samples_;
  mutable std::map<std::array<double, 5>, types::RoadSurfaceSample> contacts_;
  struct FastFootprint
  {
    std::array<double, 3> key{};
    types::RoadSurfaceSample value{};
    bool valid{false};
  };
  mutable std::array<FastFootprint, 32> fast_footprints_{};
  mutable std::size_t footprint_cursor_{0};
  mutable bool prepared_{false};
  mutable std::array<double, 3> pose_{};

public:
  explicit StageRoadCache(std::shared_ptr<const types::RoadSurface> s) : source_(std::move(s)) {}
  void prepare(double x, double y, double psi) const override
  {
    footprint_cursor_ = 0;
    std::array<double, 3> p{x, y, psi};
    if (!prepared_ || p != pose_) {
      source_->prepare(x, y, psi);
      pose_ = p;
      prepared_ = true;
      samples_.clear();
      contacts_.clear();
      for (auto & entry : fast_footprints_) entry.valid = false;
    }
  }
  types::RoadSurfaceSample sample(double x, double y) const override
  {
    auto key = std::array<double, 2>{x, y};
    auto it = samples_.find(key);
    if (it != samples_.end()) return it->second;
    return samples_.emplace(key, source_->sample(x, y)).first->second;
  }
  types::RoadSurfaceSample footprint(double x, double y, double r) const override
  {
    auto key = std::array<double, 3>{x, y, r};
    auto & fast = fast_footprints_[footprint_cursor_++ % fast_footprints_.size()];
    if (fast.valid && fast.key == key) return fast.value;
    for (const auto & entry : fast_footprints_)
      if (entry.valid && entry.key == key) {
        fast = entry;
        return fast.value;
      }
    fast = {key, source_->footprint(x, y, r), true};
    return fast.value;
  }
  types::RoadSurfaceSample tire_contact(
    double x, double y, double h, double r, double m) const override
  {
    auto key = std::array<double, 5>{x, y, h, r, m};
    auto it = contacts_.find(key);
    if (it != contacts_.end()) return it->second;
    return contacts_.emplace(key, source_->tire_contact(x, y, h, r, m)).first->second;
  }
};
struct ScopedStageRoadCache
{
  std::shared_ptr<const types::RoadSurface> & slot;
  std::shared_ptr<const types::RoadSurface> original;
  ScopedStageRoadCache(std::shared_ptr<const types::RoadSurface> & s, bool enabled)
  : slot(s), original(s)
  {
    if (enabled && s) slot = std::make_shared<StageRoadCache>(s);
  }
  ~ScopedStageRoadCache() { slot = original; }
};
}  // namespace tam::ocd::vehicle_dynamics
