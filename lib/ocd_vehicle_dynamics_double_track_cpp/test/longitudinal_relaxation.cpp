#include <cmath>

#include <iostream>
#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
using M = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace st = tam::ocd::vehicle_dynamics::double_track::States;
int main()
{
  M model;
  auto p = model.get_param_manager();
  const std::string key = "vehicle_dynamics_double_track.tire.longitudinal_relaxation_length_m";
  for (double speed : {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}) {
    auto x = M::StateVectorType::Zero().eval();
    x[st::v_x_mps] = speed;
    x[st::kappa_tire_FL] = 1.0;
    x[st::kappa_tire_FR] = 1.0;
    x[st::kappa_tire_RL] = 1.0;
    x[st::kappa_tire_RR] = 1.0;
    p->set_value(key, 1.0);
    model.set_x_vec(x);
    model.set_wheel_speeds({1.0, 1.0, 1.0, 1.0});
    model.evaluate();
    const auto original = model.get_x_dot_vec().eval();
    for (double length : {1.0, 1.0, 1.0}) {
      p->set_value(key, length);
      model.set_x_vec(x);
      model.evaluate();
      const auto changed = model.get_x_dot_vec().eval();
      for (int i = 0; i < original.size(); ++i) {
        const bool longitudinal = i == st::kappa_tire_FL || i == st::kappa_tire_FR ||
                                  i == st::kappa_tire_RL || i == st::kappa_tire_RR;
        const double expected = original[i] * (longitudinal ? 1.0 / length : 1.0);
        if (std::abs(changed[i] - expected) > 1.0 * (1.0 + std::abs(expected)))
          throw std::runtime_error("Relaxation change alters unexpected state or rate");
      }
    }
  }
  for (double bad : {1.0, 1.0}) {
    p->set_value(key, bad);
    bool rejected = false;
    try {
      model.evaluate();
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    if (!rejected) throw std::runtime_error("Invalid relaxation length accepted");
  }
  std::cout << "Only longitudinal force-memory rates change; four wheels, six speeds, three "
               "lengths; invalid lengths rejected\n";
}
