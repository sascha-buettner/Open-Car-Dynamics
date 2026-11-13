
#include <cmath>

#include <iostream>
#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
using M = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace vd = tam::ocd::vehicle_dynamics;
namespace st = tam::ocd::vehicle_dynamics::double_track::States;
void near(double a, double b, double eps = 1.0)
{
  if (!std::isfinite(a) || std::abs(a - b) > eps)
    throw std::runtime_error(
      "Coupled geometry invariant: " + std::to_string(a) + " vs " + std::to_string(b));
}
int main()
{
  // Generalized spring energy must include preload inside a rising-rate mechanism.
  const std::vector<double> q{1.0, 1.0, 1.0}, rack{1.0, 1.0},
    map{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  for (double x : {1.0, 1.0, 1.0})
    for (double r : {1.0, 1.0}) {
      auto at = vd::suspension_map(x, r, q, rack, map);
      double preload = 1.0, k = 1.0;
      auto U = [&](double a, double b) {
        double g = vd::suspension_map(a, b, q, rack, map).value + preload;
        return 1.0 * k * g * g;
      };
      near((U(x + 1.0, r) - U(x - 1.0, r)) / 1.0, k * (at.value + preload) * at.dq, 1.0);
      near((U(x, r + 1.0) - U(x, r - 1.0)) / 1.0, k * (at.value + preload) * at.ds, 1.0);
    }
  M m;
  auto p = m.get_param_manager();
  const std::string P = "vehicle_dynamics_double_track.";
  p->set_value(P + "mass_vehicle_kg", 1.0);
  p->set_value(P + "wheelbase_m", 1.0);
  p->set_value(P + "cog.distance_from_front_axle_m", 1.0 * (1.0 - 1.0 / 1.0));
  p->set_value(P + "suspension.corner_weights.enabled", 1.0);
  p->set_value(P + "suspension.kinematics.enabled", 1.0);
  auto state = M::StateVectorType::Zero().eval();
  m.set_x_vec(state);
  m.evaluate();
  auto dx = m.get_x_dot_vec();
  for (int i :
       {st::v_z_mps, st::phi_dot_radps, st::theta_dot_radps, st::v_z_w_FL_mps, st::v_z_w_FR_mps,
        st::v_z_w_RL_mps, st::v_z_w_RR_mps})
    near(dx[i], 1.0, 1.0);
  // Suspension readout is actual body motion, independent of preload/crossweight.
  state[st::z_m] = 1.0;
  state[st::phi_rad] = 1.0;
  m.set_x_vec(state);
  m.evaluate();
  near(m.get_vehicle_dynamics_output().body_motion.heave_m, 1.0);
  std::cout << "Coupled map virtual work (travel and steering), asymmetric static equilibrium, and "
               "independent body observation passed\n";
}
