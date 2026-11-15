#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
#include "synthetic_vehicle.hpp"

using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace st = tam::ocd::vehicle_dynamics::double_track::States;

double scalar(const tam::tsl::data_map_t & data, const std::string & key)
{
  return std::visit(
    [](const auto & v) -> double {
      using T = std::decay_t<decltype(v)>;
      if constexpr (std::is_arithmetic_v<T>) return v;
      else if constexpr (std::is_pointer_v<T>) return *v;
      else throw std::runtime_error("Non-scalar ring log");
    }, data.at(key));
}

int main()
{
  Model model;
  synthetic_vehicle::contact(model, true);
  auto pm = model.get_param_manager();
  const std::string base = "vehicle_dynamics_double_track.tire.mfschwift.";
  for (const auto & item : {
         std::pair<const char *, double>{"inertia_transverse_kgm2", .04},
         {"inertia_spin_kgm2", .08}, {"frequency_lateral_hz", 40.},
         {"frequency_yaw_camber_hz", 35.}, {"frequency_windup_hz", 45.},
         {"damping_lateral", .15}, {"damping_yaw_camber", .15},
         {"damping_windup", .15}, {"longitudinal_stiffness_Npm", 100000.}})
    pm->set_value(base + item.first, item.second);
  pm->set_value(base + "enabled", 1.);
  auto state = Model::StateVectorType::Zero().eval();
  model.set_x_vec(state);
  model.evaluate();
  if (!model.get_x_dot_vec().allFinite())
    throw std::runtime_error("Synthetic ring equilibrium is nonfinite");
  state[st::ring_x_FL_m] = .001;
  model.set_x_vec(state);
  model.evaluate();
  const double spring_accel = model.get_x_dot_vec()[st::ring_vx_FL_mps];
  const double energy = scalar(model.get_logger()->get_data(), "mfschwift/ring_energy_J/front_left");
  if (!(spring_accel < 0. && energy > 0.))
    throw std::runtime_error("Ring sidewall does not restore displacement");
  state[st::ring_vx_FL_mps] = .2;
  model.set_x_vec(state);
  model.evaluate();
  const double damping = scalar(model.get_logger()->get_data(), "mfschwift/sidewall_dissipation_W/front_left");
  if (!(std::isfinite(damping) && damping > 0.))
    throw std::runtime_error("Ring sidewall is not dissipative");
  pm->set_value(base + "inertia_spin_kgm2", 0.);
  bool rejected = false;
  try { model.evaluate(); } catch (const std::invalid_argument &) { rejected = true; }
  if (!rejected) throw std::runtime_error("Zero ring inertia accepted");
}
