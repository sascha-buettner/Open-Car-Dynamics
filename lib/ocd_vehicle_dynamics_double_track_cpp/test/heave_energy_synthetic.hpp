#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"

namespace heave_energy_synthetic
{
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace st = tam::ocd::vehicle_dynamics::double_track::States;

double value(const tam::tsl::data_map_t & data, const std::string & key)
{
  return std::visit(
    [](const auto & v) -> double {
      using T = std::decay_t<decltype(v)>;
      if constexpr (std::is_arithmetic_v<T>) return v;
      else if constexpr (std::is_pointer_v<T>) return *v;
      else throw std::runtime_error("Non-scalar heave log");
    }, data.at(key));
}

void run(bool coil)
{
  Model model;
  auto parameters = model.get_param_manager();
  const std::string base = "vehicle_dynamics_double_track.suspension.heave.";
  const std::string stiffness = coil ? "spring_stiffness_Npm." : "rubber_stiffness_Npm2.";
  parameters->set_value(base + stiffness + "front", coil ? 40000. : 1000000.);
  parameters->set_value(base + stiffness + "rear", coil ? 30000. : 800000.);
  parameters->set_value(base + "device_gap_m.front", 0.);
  parameters->set_value(base + "device_gap_m.rear", 0.);
  auto state = Model::StateVectorType::Zero().eval();
  const auto evaluate = [&](const auto & s) {
    model.set_x_vec(s);
    model.evaluate();
    return model.get_logger()->get_data();
  };
  const auto initial = evaluate(state);
  for (const char * wheel : {"front_left", "rear_left"})
    if (std::abs(value(initial, std::string("imr/suspension_heave_force_N/") + wheel)) > 1e-8)
      throw std::runtime_error("Third element changed static preload");
  state[st::z_w_FL_m] = .01;
  state[st::z_w_FR_m] = .01;
  state[st::z_w_RL_m] = .015;
  state[st::z_w_RR_m] = .015;
  const auto data = evaluate(state);
  for (const auto & item : {std::pair<int, const char *>{st::z_w_FL_m, "front_left"},
                           {st::z_w_RR_m, "rear_right"}}) {
    const double force = value(data, std::string("imr/suspension_heave_force_N/") + item.second);
    const auto total_energy = [&](const auto & s) {
      const auto d = evaluate(s);
      double sum = 0.;
      for (const char * wheel : {"front_left", "front_right", "rear_left", "rear_right"})
        sum += value(d, std::string("imr/suspension_heave_energy_J/") + wheel);
      return sum;
    };
    const double step = 1e-6;
    auto plus = state, minus = state;
    plus[item.first] += step;
    minus[item.first] -= step;
    const double gradient = (total_energy(plus) - total_energy(minus)) / (2. * step);
    if (!(std::isfinite(force) && std::isfinite(gradient) &&
          std::abs(force + gradient) < .001 * std::max(1., std::abs(gradient))))
      throw std::runtime_error("Third-element wheel force violates virtual work");
  }
  auto roll = Model::StateVectorType::Zero().eval();
  roll[st::z_w_FL_m] = .01;
  roll[st::z_w_FR_m] = -.01;
  const auto rolled = evaluate(roll);
  if (std::abs(value(rolled, "imr/suspension_heave_force_N/front_left")) > 1e-8)
    throw std::runtime_error("Pure axle roll excites third element");
}
}  // namespace heave_energy_synthetic
