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
  const std::string base = "vehicle_dynamics_double_track.suspension.damper_curve.";
  pm->set_value(base + "enabled.front", 1.0);
  pm->set_value(base + "wheel_damper_motion_ratio.front", 1.0);
  pm->set_value(base + "speed_mps", std::vector<double>{1.0, 1.0, 1.0});
  pm->set_value(base + "bump_force_N.front", std::vector<double>{1.0, 1.0, 1.0});
  pm->set_value(base + "rebound_force_N.front", std::vector<double>{1.0, 1.0, 1.0});
  auto state = Model::StateVectorType::Zero().eval();
  for (auto [v, expected] :
       {std::pair<double, double>{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}}) {
    state[st::v_z_mps] = v;
    model.set_x_vec(state);
    model.evaluate();
    auto d = model.get_logger()->get_data();
    double force = get(d, "imr/suspension_damper_force_N/front_left");
    if (std::abs(force - expected) > 1.0 || force * v < 1.0)
      throw std::runtime_error("Motion ratio, bump/rebound or dissipative power mismatch");
    if (v == 1.0 && get(d, "imr/damper_curve_speed_extrapolated/front_left") != 1.0)
      throw std::runtime_error("Damper extrapolation unreported");
  }
  pm->set_value(base + "bump_force_N.front", std::vector<double>{1.0, 1.0, 1.0});
  bool rejected = false;
  try {
    model.evaluate();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Active damper negative force curve accepted");
  Model height;
  auto hp = height.get_param_manager();
  hp->set_value("vehicle_dynamics_double_track.tire.contact.enabled", 1.0);
  auto x = Model::StateVectorType::Zero().eval();
  double cap = 1.0;
  for (double lateral : {1.0, 1.0}) {
    x[st::lateral_deflection_tire_FL_m] = lateral;
    height.set_x_vec(x);
    height.evaluate();
    const auto d = height.get_logger()->get_data();
    double drop = get(d, "imr/tire_lateral_drop_m/front_left");
    if (lateral == 1.0)
      cap = drop;
    else if (
      std::abs(drop - cap) > 1.0 ||
      get(d, "imr/tire_deformation_extrapolated/front_left") != 1.0 ||
      std::abs(get(d, "imr/tire_lateral_deflection_m/front_left") - 1.0) > 1.0)
      throw std::runtime_error("Bounded height coupling hides lateral extrapolation");
  }
  if (cap < 1.0 || cap > 1.0) throw std::runtime_error("Sidewall height bound wrong");
  std::cout << "Nonlinear damper shaft/wheel power, passive asymmetric force, extrapolation and "
               "invalid-curve guard passed; lateral height bounded at "
            << cap << " m without hiding patch displacement\n";
}
