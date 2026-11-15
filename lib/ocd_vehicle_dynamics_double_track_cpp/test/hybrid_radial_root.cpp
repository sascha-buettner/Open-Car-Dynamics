#include <stdexcept>
#include "ocd_vehicle_dynamics_double_track_cpp/hybrid_radial_root.hpp"
using tam::ocd::vehicle_dynamics::hybrid_radial_root;
int main()
{
  const Eigen::Vector3d exact(2., 3., 4.), lower(1., 1., 1.), upper(5., 5., 5.);
  Eigen::Matrix3d a;
  a << 2., .2, .1, .1, 3., .2, .2, .1, 4.;
  const auto root = hybrid_radial_root(
    [&](const Eigen::Vector3d &x) { return (a * (x - exact)).eval(); },
    Eigen::Vector3d(1.5, 2., 3.), lower, upper);
  if ((root - exact).norm() > 1e-7) throw std::runtime_error("Hybrid coupled root incorrect");
  for (int kind : {0, 1}) {
    bool rejected = false;
    try {
      (void)hybrid_radial_root(
        [&](const Eigen::Vector3d &x) {
          return Eigen::Vector3d(kind == 0 ? x[0] + 1. : (x[0] < 2. ? -1. : 1.),
                                 x[1] - 3., x[2] - 4.);
        }, Eigen::Vector3d(1.5, 2., 3.), lower, upper);
    } catch (const std::runtime_error &) { rejected = true; }
    if (!rejected) throw std::runtime_error("Hybrid accepted missing or discontinuous root");
  }
}
