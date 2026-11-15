#include <cmath>
#include <stdexcept>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
#include "synthetic_vehicle.hpp"
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace states = tam::ocd::vehicle_dynamics::double_track::States;

double scalar(const tam::tsl::data_map_t & data, const std::string & key)
{
  return std::visit(
    [](const auto & v) -> double {
      using T = std::decay_t<decltype(v)>;
      if constexpr (std::is_arithmetic_v<T>) return v;
      else if constexpr (std::is_pointer_v<T>) return *v;
      else return NAN;
    }, data.at(key));
}

int main()
{
  Model model;
  synthetic_vehicle::contact(model, true);
  auto state = Model::StateVectorType::Zero().eval();
  model.set_x_vec(state);
  model.evaluate();
  if (!model.get_x_dot_vec().allFinite())
    throw std::runtime_error("Synthetic belt equilibrium is not finite");
  const std::string force_key = "imr/belt_carcass_force_N/front_left";
  const std::string energy_key = "imr/belt_elastic_energy_J/front_left";
  const double initial = scalar(model.get_logger()->get_data(), force_key);
  state[states::z_belt_FL_m] = .001;
  model.set_x_vec(state);
  model.evaluate();
  const double displaced = scalar(model.get_logger()->get_data(), force_key);
  const double energy = scalar(model.get_logger()->get_data(), energy_key);
  if (!(displaced > initial && energy >= 0.))
    throw std::runtime_error("Belt elastic force or energy did not respond to displacement");
  state[states::v_z_belt_FL_mps] = .1;
  model.set_x_vec(state);
  model.evaluate();
  const double damped = scalar(model.get_logger()->get_data(), force_key);
  if (!(damped > displaced && (damped - displaced) * .1 > 0.))
    throw std::runtime_error("Belt damping is not passive");
  model.get_param_manager()->set_value(
    "vehicle_dynamics_double_track.tire.contact.belt.mass_kg.front", 15.);
  bool rejected = false;
  try {
    model.evaluate();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Belt mass equal to unsprung mass accepted");
}
