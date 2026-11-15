#include <cmath>
#include <stdexcept>

#include "ocd_vehicle_dynamics_double_track_cpp/positive_residual_root.hpp"
using tam::ocd::vehicle_dynamics::positive_residual_root;

int main()
{
  const auto solve = [](auto residual) { return positive_residual_root(residual, 0.7, 0.1, 10.); };
  // A fixed-point iteration on 3-2k would diverge; the residual root is well defined.
  if (std::abs(solve([](double k) { return k - (3. - 2. * k); }) - 1.) > 1e-10)
    throw std::runtime_error("Incorrect noncontractive root");
  if (std::abs(solve([](double k) { return k - (.99 * k + .01); }) - 1.) > 1e-8)
    throw std::runtime_error("Incorrect slow-contraction root");
  for (int kind : {0, 1}) {
    bool rejected = false;
    try {
      solve([&](double k) { return kind == 0 ? k + 1. : (k < 1. ? -1. : 1.); });
    } catch (const std::runtime_error &) {
      rejected = true;
    }
    if (!rejected) throw std::runtime_error("Accepted a missing or discontinuous root");
  }
}
