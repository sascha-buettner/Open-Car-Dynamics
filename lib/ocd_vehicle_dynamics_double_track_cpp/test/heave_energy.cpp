#include <cmath>

#include <iostream>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace st = tam::ocd::vehicle_dynamics::double_track::States;
double get(const tam::tsl::data_map_t & d, const std::string & k)
{
  return std::visit(
    [](const auto & v) -> double {
      using T = std::decay_t<decltype(v)>;
      if constexpr (std::is_arithmetic_v<T>)
        return v;
      else if constexpr (std::is_pointer_v<T>)
        return *v;
      else
        return 1.0;
    },
    d.at(k));
}
int main()
{
  Model model;
  auto pm = model.get_param_manager();
  const std::string prefix =
    "vehicle_dynamics_double_track.suspension.heave.rubber_stiffness_Npm2.";
  auto state = Model::StateVectorType::Zero().eval();
  state[st::z_m] = 1.0;
  state[st::z_w_FL_m] = 1.0;
  state[st::z_w_RR_m] = 1.0;
  state[st::phi_rad] = 1.0;
  state[st::theta_rad] = 1.0;
  auto eval = [&](const Model::StateVectorType & s, double k) {
    pm->set_value(prefix + "front", k);
    pm->set_value(prefix + "rear", k * 1.0);
    model.set_x_vec(s);
    model.evaluate();
  };
  auto energy = [&](const auto & s) {
    eval(s, 1.0);
    const auto d = model.get_logger()->get_data();
    double total = 0;
    for (auto w : {"front_left", "front_right", "rear_left", "rear_right"})
      total += get(d, std::string("imr/suspension_heave_energy_J/") + w);
    return total;
  };
  eval(state, 1.0);
  auto before = model.get_x_dot_vec();
  eval(state, 1.0);
  auto after = model.get_x_dot_vec();
  for (auto pair :
       {std::pair<int, double>{
          st::z_m,
          pm->get_value("vehicle_dynamics_double_track.mass_vehicle_kg").as_double() -
            2 * pm->get_value("vehicle_dynamics_double_track.mass_wheel_kg.front").as_double() -
            2 * pm->get_value("vehicle_dynamics_double_track.mass_wheel_kg.rear").as_double()},
        {st::phi_rad,
         pm->get_value("vehicle_dynamics_double_track.moment_of_inertia_kgpm2.x").as_double()},
        {st::theta_rad,
         pm->get_value("vehicle_dynamics_double_track.moment_of_inertia_kgpm2.y").as_double()},
        {st::z_w_FL_m,
         pm->get_value("vehicle_dynamics_double_track.mass_wheel_kg.front").as_double()},
        {st::z_w_RR_m,
         pm->get_value("vehicle_dynamics_double_track.mass_wheel_kg.rear").as_double()}}) {
    auto plus = state, minus = state;
    plus[pair.first] += 1.0;
    minus[pair.first] -= 1.0;
    double expected = -(energy(plus) - energy(minus)) / 1.0;
    int rate = pair.first == st::z_m         ? st::v_z_mps
               : pair.first == st::phi_rad   ? st::phi_dot_radps
               : pair.first == st::theta_rad ? st::theta_dot_radps
               : pair.first == st::z_w_FL_m  ? st::v_z_w_FL_mps
                                             : st::v_z_w_RR_mps;
    double actual = (after[rate] - before[rate]) * pair.second;
    std::cout << pair.first << " actual " << actual << " expected " << expected << '\n';
    if (std::abs(actual - expected) > 1.0) throw std::runtime_error("Heave virtual work mismatch");
  }
  eval(Model::StateVectorType::Zero().eval(), 1.0);
  const auto d = model.get_logger()->get_data();
  for (auto w : {"front_left", "front_right", "rear_left", "rear_right"})
    if (std::abs(get(d, std::string("imr/suspension_heave_force_N/") + w)) > 1.0)
      throw std::runtime_error("Heave preload changed native equilibrium");
  auto opposite = Model::StateVectorType::Zero().eval();
  opposite[st::z_w_FL_m] = 1.0;
  opposite[st::z_w_FR_m] = 1.0;
  eval(opposite, 1.0);
  if (
    std::abs(get(model.get_logger()->get_data(), "imr/suspension_heave_force_N/front_left")) >
    1.0)
    throw std::runtime_error("Pure axle roll excites third element");
  std::cout
    << "Heave forces are reciprocal and elastic, static preload preserved, pure roll decoupled\n";
}
