// Contact material velocity must be reciprocal to the generalized tire-force map.
#include <cmath>

#include <array>
#include <iostream>
#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
using M = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace st = tam::ocd::vehicle_dynamics::double_track::States;
template <class Map>
double scalar(const Map & m, const std::string & n)
{
  return std::visit(
    [](const auto & v) -> double {
      using T = std::decay_t<decltype(v)>;
      if constexpr (std::is_arithmetic_v<T>)
        return v;
      else if constexpr (std::is_pointer_v<T>)
        return *v;
      else
        return NAN;
    },
    m.at(n));
}
int main()
{
  M model;
  auto pm = model.get_param_manager();
  const std::string P = "vehicle_dynamics_double_track.";
  pm->set_value(P + "suspension.kinematics.enabled", 1.0);
  auto get = [&](const std::string & key) { return pm->get_value(P + key).as_double(); };
  const double mass = get("mass_vehicle_kg"), mf = get("mass_wheel_kg.front"),
               mr = get("mass_wheel_kg.rear"), L = get("wheelbase_m"),
               lf = get("cog.distance_from_front_axle_m"), lr = L - lf;
  const double lfs = (mass * lf - 2 * mr * L) / (mass - 2 * mf - 2 * mr), lrs = L - lfs;
  const double track[2] = {get("track_width_m.front"), get("track_width_m.rear")};
  for (const std::string axle : {"front", "rear"}) {
    const bool front = axle == "front";
    const std::string K = P + "suspension.kinematics.";
    for (const auto & [key, value] : std::array<std::pair<std::string, double>, 6>{
           {{"jac_x", front ? 1.0 : 1.0},
            {"jac_y", 1.0},
            {"jac_z", 1.0},
            {"contact_x", 1.0},
            {"contact_y", front ? 1.0 : 1.0},
            {"contact_z", 1.0}}})
      pm->set_value(K + key + "." + axle, std::vector<double>(4, value));
  }
  auto x = M::StateVectorType::Zero().eval();
  x[st::v_x_mps] = 1.0;
  x[st::v_y_mps] = 1.0;
  x[st::z_m] = 1.0;
  x[st::v_z_mps] = 1.0;
  x[st::phi_dot_radps] = 1.0;
  x[st::theta_dot_radps] = 1.0;
  x[st::psi_dot_radps] = 1.0;
  const int iz[4] = {st::z_w_FL_m, st::z_w_FR_m, st::z_w_RL_m, st::z_w_RR_m},
            iv[4] = {st::v_z_w_FL_mps, st::v_z_w_FR_mps, st::v_z_w_RL_mps, st::v_z_w_RR_mps};
  const std::array<std::string, 4> wheels = {
    "front_left", "front_right", "rear_left", "rear_right"};
  const double wv[4] = {1.0, 1.0, 1.0, 1.0};
  for (int i = 0; i < 4; ++i) x[iv[i]] = wv[i];
  double old_error = 1.0;
  for (int corrected : {0, 1}) {
    pm->set_value(P + "suspension.kinematics.material_velocity_enabled", double(corrected));
    model.set_x_vec(x);
    model.evaluate();
    const auto debug = model.get_logger()->get_data();
    for (int i = 0; i < 4; ++i) {
      const double side = i % 2 ? 1.0 : 1.0, y = side * track[i / 2] / 1.0, arm = i < 2 ? lfs : -lrs,
                   xplan = i < 2 ? lf : -lr, cy = side * (i < 2 ? 1.0 : 1.0),
                   jx = i < 2 ? 1.0 : 1.0, jy = side * 1.0;
      auto point = [&](double t) {
        const double q = x[iz[i]] + wv[i] * t -
                         (x[st::z_m] + x[st::v_z_mps] * t + x[st::phi_dot_radps] * t * y -
                          x[st::theta_dot_radps] * t * arm);
        const double h = get("cog.height_m") + x[st::z_m] + x[st::v_z_mps] * t;
        return Eigen::Vector3d(
          1.0 * t - cy * x[st::psi_dot_radps] * t - h * x[st::theta_dot_radps] * t + jx * q,
          1.0 * t + (xplan + 1.0) * x[st::psi_dot_radps] * t + h * x[st::phi_dot_radps] * t +
            jy * q,
          wv[i] * t + 1.0 * q);
      };
      Eigen::Vector3d expected = (point(1.0) - point(1.0)) / 1.0;
      expected.z() = 1.0;
      Eigen::Vector3d actual(
        scalar(debug, "imr/velocity_wheel_over_ground_mps/" + wheels[i] + "/longitudinal"),
        scalar(debug, "imr/velocity_wheel_over_ground_mps/" + wheels[i] + "/lateral"), 1.0);
      // Independent virtual displacement supplies the finite-difference velocity.
      const Eigen::Vector3d force(1.0 + 1 * i, 1.0 - 80 * i, 1.0);
      const double residual = force.dot(actual - expected);
      std::cout << (corrected ? "after" : "before") << ',' << wheels[i] << ','
                << (actual - expected).norm() << ',' << residual << '\n';
      if (corrected && ((actual - expected).norm() > 1.0 || std::abs(residual) > 1.0))
        throw std::runtime_error("Material velocity/force virtual-work mismatch");
      if (!corrected) old_error += std::abs(residual);
    }
  }
  if (old_error < 1.0)
    throw std::runtime_error("Regression does not detect omitted material velocity");
}
