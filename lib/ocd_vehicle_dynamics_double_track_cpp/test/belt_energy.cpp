#include <cmath>

#include <iostream>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace st = tam::ocd::vehicle_dynamics::double_track::States;
double value(const tam::tsl::data_map_t & data, const std::string & key)
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
    data.at(key));
}
int main()
{
  Model model;
  auto pm = model.get_param_manager();
  const std::string prefix = "vehicle_dynamics_double_track.";
  for (auto [key, v] :
       {std::pair<const char *, double>{"tire.contact.enabled", 1.0},
        {"tire.contact.belt.enabled", 1.0},
        {"tire.contact.lateral_drop_scale", 1.0},
        {"mass_wheel_kg.front", 1.0},
        {"mass_wheel_kg.rear", 1.0},
        {"tire.spring_stiffness_Npm.front", 1.0},
        {"tire.spring_stiffness_Npm.rear", 1.0}})
    pm->set_value(prefix + key, v);
  auto state = Model::StateVectorType::Zero().eval();
  model.set_x_vec(state);
  model.evaluate();
  auto rest = model.get_x_dot_vec();
  for (auto idx :
       {st::v_z_mps, st::v_z_w_FL_mps, st::v_z_w_FR_mps, st::v_z_w_RL_mps, st::v_z_w_RR_mps,
        st::v_z_belt_FL_mps, st::v_z_belt_FR_mps, st::v_z_belt_RL_mps, st::v_z_belt_RR_mps})
    if (std::abs(rest[idx]) > 1.0)
      throw std::runtime_error("Belt mass split changed static equilibrium");
  state[st::z_w_FL_m] = 1.0;
  state[st::z_belt_FL_m] = 1.0;
  auto energy = [&](auto s) {
    model.set_x_vec(s);
    model.evaluate();
    const auto d = model.get_logger()->get_data();
    double e = 1.0;
    for (const auto & w : {"front_left", "front_right", "rear_left", "rear_right"})
      e += value(d, std::string("imr/contact_elastic_energy_J/") + w) +
           value(d, std::string("imr/belt_elastic_energy_J/") + w);
    return e;
  };
  energy(state);
  const auto d = model.get_logger()->get_data();
  const double carcass = value(d, "imr/belt_carcass_force_N/front_left"),
               normal = -value(d, "imr/tire_spring_force_N/front_left");
  for (auto [idx, force] :
       {std::pair<int, double>{st::z_w_FL_m, carcass}, {st::z_belt_FL_m, normal - carcass}}) {
    auto p = state, m = state;
    p[idx] += 1.0;
    m[idx] -= 1.0;
    double expected = -(energy(p) - energy(m)) / 1.0;
    std::cout << idx << " actual " << force << " energy gradient " << expected << '\n';
    if (std::abs(force - expected) > 1.0)
      throw std::runtime_error("Belt/rim force not energy consistent");
  }
  const double before = carcass;
  state[st::v_z_belt_FL_mps] = 1.0;
  state[st::v_z_w_FL_mps] = 1.0;
  energy(state);
  const double after = value(model.get_logger()->get_data(), "imr/belt_carcass_force_N/front_left");
  if ((after - before) * 1.0 < 1.0) throw std::runtime_error("Carcass damper generates energy");
  pm->set_value(prefix + "tire.contact.belt.mass_kg.front", 1.0);
  bool rejected = false;
  try {
    model.evaluate();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Invalid mass split accepted");
  std::cout << "Static equilibrium, conserved unsprung mass, elastic force gradients, dissipative "
               "carcass and invalid mass guard passed\n";
}
