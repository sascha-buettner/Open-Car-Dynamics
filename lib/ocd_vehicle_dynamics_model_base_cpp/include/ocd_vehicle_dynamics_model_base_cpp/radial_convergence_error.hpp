#pragma once

#include <stdexcept>

namespace tam::ocd::vehicle_dynamics
{
class RadialConvergenceError : public std::runtime_error
{
public:
  using std::runtime_error::runtime_error;
};
}  // namespace tam::ocd::vehicle_dynamics
