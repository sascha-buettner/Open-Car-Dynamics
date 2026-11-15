#include <cmath>
#include <stdexcept>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/suspension_lookup.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
#include "synthetic_vehicle.hpp"
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace states = tam::ocd::vehicle_dynamics::double_track::States;

int main()
{
  using tam::ocd::vehicle_dynamics::suspension_lookup;
  const auto sample = suspension_lookup(.05, {-.1, 0., .1}, {-.02, 0., .02});
  if (std::abs(sample.first - .01) > 1e-12 || std::abs(sample.second - .2) > 1e-12)
    throw std::runtime_error("Suspension lookup interpolation failed");
  bool rejected = false;
  try {
    suspension_lookup(.05, {0., 0., .1}, {0., 0., .02});
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Unsorted suspension lookup accepted");

  Model model;
  auto state = Model::StateVectorType::Zero().eval();
  model.set_x_vec(state);
  model.evaluate();
  const auto before = model.get_x_dot_vec().eval();
  synthetic_vehicle::geometry(model);
  model.evaluate();
  if (!model.get_x_dot_vec().allFinite() ||
      std::abs(model.get_x_dot_vec()[states::v_z_mps] - before[states::v_z_mps]) > 1e-8)
    throw std::runtime_error("Flat geometry changed static body equilibrium");
  const auto parameters = model.get_param_manager();
  parameters->set_value("vehicle_dynamics_double_track.suspension.geometry.toe_delta_rad.front",
                        std::vector<double>{-.01, 0., .01});
  state[states::z_w_FL_m] = .05;
  state[states::z_m] = .01;
  model.set_x_vec(state);
  model.evaluate();
  const auto steering = model.get_vehicle_dynamics_output();
  if (!model.get_x_dot_vec().allFinite() ||
      std::abs(steering.body_motion.heave_m - .01) > 1e-12)
    throw std::runtime_error("Geometry contaminated body heave observation");
  const auto data = model.get_logger()->get_data();
  const auto it = data.find("imr/effective_steering_angle_per_wheel_rad/front_left");
  if (it == data.end()) throw std::runtime_error("Steering geometry was not logged");
  const double toe = std::visit(
    [](const auto & v) -> double {
      using T = std::decay_t<decltype(v)>;
      if constexpr (std::is_arithmetic_v<T>) return v;
      else if constexpr (std::is_pointer_v<T>) return *v;
      else return NAN;
    }, it->second);
  if (!(std::isfinite(toe) && std::abs(toe) > 1e-4))
    throw std::runtime_error("Configured toe geometry had no effect");
}
