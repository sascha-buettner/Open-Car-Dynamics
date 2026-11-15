#include <cmath>
#include <stdexcept>
#include "ocd_vehicle_dynamics_double_track_cpp/coupled_radial_root.hpp"
using tam::ocd::vehicle_dynamics::coupled_radial_root;
int main()
{
  const Eigen::Vector3d exact(2., 3., 4.), lower(1., 1., 1.), upper(5., 5., 5.);
  Eigen::Matrix3d a;
  a << 2., .2, .1, .1, 3., .2, .2, .1, 4.;
  const auto linear = coupled_radial_root(
    [&](const Eigen::Vector3d &x) { return (a * (x - exact)).eval(); },
    Eigen::Vector3d(1.5, 2., 3.), lower, upper);
  if ((linear - exact).norm() > 1e-7) throw std::runtime_error("Linear coupled root incorrect");
  const auto nonlinear = coupled_radial_root(
    [&](const Eigen::Vector3d &x) {
      return Eigen::Vector3d(x[0] * x[0] - 4. + .1 * (x[1] - 3.),
                             x[1] * x[1] - 9. + .1 * (x[2] - 4.),
                             x[2] * x[2] - 16. + .1 * (x[0] - 2.));
    }, Eigen::Vector3d(1.8, 2.8, 3.8), lower, upper);
  if ((nonlinear - exact).norm() > 1e-7) throw std::runtime_error("Nonlinear coupled root incorrect");
  for (int kind : {0, 1}) {
    bool rejected = false;
    try {
      (void)coupled_radial_root(
        [&](const Eigen::Vector3d &x) {
          return Eigen::Vector3d(kind == 0 ? x[0] + 1. : (x[0] < 2. ? -1. : 1.),
                                 x[1] - 3., x[2] - 4.);
        }, Eigen::Vector3d(1.5, 2., 3.), lower, upper);
    } catch (const std::runtime_error &) { rejected = true; }
    if (!rejected) throw std::runtime_error("Accepted missing or discontinuous root");
  }
}
