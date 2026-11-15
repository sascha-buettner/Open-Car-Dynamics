#include <cmath>
#include <stdexcept>
#include <vector>
#include "ocd_vehicle_dynamics_double_track_cpp/coupled_radial_root.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/generic_radial_compliance.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/hybrid_radial_root.hpp"

using tam::ocd::vehicle_dynamics::GenericRadialCompliance;
using tam::ocd::vehicle_dynamics::GenericRadialLookup;
using tam::ocd::vehicle_dynamics::coupled_radial_root;
using tam::ocd::vehicle_dynamics::hybrid_radial_root;

int main()
{
  // Synthetic grid and coupled deflections. No measured tire values are embedded.
  GenericRadialLookup lookup;
  std::vector<double> front(96, 10000.), rear(96, 12000.);
  lookup.configure({1., 2., 3.}, {1., 2., 3., 4.}, {1., 2., 3., 4., 5., 6., 7., 8.}, front, rear);
  GenericRadialCompliance compliance;
  compliance.reset(lookup, false, 2., 2.);
  const Eigen::Vector3d exact(2., 3., 4.), lower(1., 1., 1.), upper(5., 5., 5.);
  const auto make_residual = [&](const Eigen::Vector3d &x) {
    Eigen::Vector3d r;
    for (int i = 0; i < 3; ++i) {
      const double displacement = compliance.values(1000. * x[i])[0];
      const double target = compliance.values(1000. * exact[i])[0];
      r[i] = displacement - target + .01 * (x[(i + 1) % 3] - exact[(i + 1) % 3]);
    }
    return r;
  };
  const auto root = coupled_radial_root(make_residual, Eigen::Vector3d(1.5, 2.5, 3.5), lower, upper);
  if ((root - exact).norm() > 1e-7 || make_residual(root).lpNorm<Eigen::Infinity>() > 1e-9)
    throw std::runtime_error("Coupled radial force/compliance residual not satisfied");
  const auto hybrid = hybrid_radial_root(make_residual, Eigen::Vector3d(1.5, 2.5, 3.5), lower, upper);
  if ((hybrid - exact).norm() > 1e-7 || make_residual(hybrid).lpNorm<Eigen::Infinity>() > 1e-9)
    throw std::runtime_error("Hybrid radial force/compliance residual not satisfied");
}
