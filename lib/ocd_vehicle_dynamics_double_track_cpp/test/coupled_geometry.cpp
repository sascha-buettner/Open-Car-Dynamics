#include <cmath>
#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/suspension_table_2d.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
#include "synthetic_vehicle.hpp"
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace vd = tam::ocd::vehicle_dynamics;
namespace states = tam::ocd::vehicle_dynamics::double_track::States;

int main()
{
  const std::vector<double> travel{-.1, .1}, rack{-.1, .1};
  const std::vector<double> map{-.11, .09, -.09, .11};
  const auto at = vd::suspension_map(.03, .04, travel, rack, map);
  if (std::abs(at.value - .034) > 1e-12 ||
      std::abs(at.dq - 1.) > 1e-12 || std::abs(at.ds - .1) > 1e-12)
    throw std::runtime_error("Coupled jounce/rack map derivative incorrect");
  bool rejected = false;
  try {
    vd::suspension_map(0., 0., travel, {0., 0.}, map);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Invalid rack knots accepted");

  Model model;
  synthetic_vehicle::kinematics(model);
  auto parameters = model.get_param_manager();
  const std::string base = "vehicle_dynamics_double_track.suspension.corner_weights.";
  parameters->set_value(base + "kg.front_left", 205.);
  parameters->set_value(base + "kg.front_right", 195.);
  parameters->set_value(base + "kg.rear_left", 195.);
  parameters->set_value(base + "kg.rear_right", 205.);
  parameters->set_value(base + "enabled", 1.);
  auto state = Model::StateVectorType::Zero().eval();
  state[states::z_m] = .01;
  model.set_x_vec(state);
  model.evaluate();
  if (!model.get_x_dot_vec().allFinite() ||
      std::abs(model.get_vehicle_dynamics_output().body_motion.heave_m - .01) > 1e-12)
    throw std::runtime_error("Synthetic coupled geometry fails finite body response");
  parameters->set_value(base + "kg.front_left", 204.);
  rejected = false;
  try {
    model.evaluate();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Nonconserving corner weights accepted");
}
