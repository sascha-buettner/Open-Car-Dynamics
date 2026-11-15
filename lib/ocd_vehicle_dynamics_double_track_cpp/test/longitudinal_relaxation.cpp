#include <cmath>
#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace states = tam::ocd::vehicle_dynamics::double_track::States;

int main()
{
  Model model;
  auto parameters = model.get_param_manager();
  const std::string key = "vehicle_dynamics_double_track.tire.longitudinal_relaxation_length_m";
  for (double speed : {0., 5., 10., 20., 50., 100.}) {
    auto state = Model::StateVectorType::Zero().eval();
    state[states::v_x_mps] = speed;
    for (int i : {states::kappa_tire_FL, states::kappa_tire_FR,
                  states::kappa_tire_RL, states::kappa_tire_RR})
      state[i] = -.1;
    model.set_x_vec(state);
    model.set_wheel_speeds({speed / .3, speed / .3, speed / .3, speed / .3});
    parameters->set_value(key, .3);
    model.evaluate();
    const auto baseline = model.get_x_dot_vec().eval();
    for (double length : {.15, .6}) {
      parameters->set_value(key, length);
      model.evaluate();
      const auto changed = model.get_x_dot_vec().eval();
      for (int i = 0; i < changed.size(); ++i) {
        const bool slip = i == states::kappa_tire_FL || i == states::kappa_tire_FR ||
                          i == states::kappa_tire_RL || i == states::kappa_tire_RR;
        const double expected = baseline[i] * (slip ? .3 / length : 1.);
        if (std::abs(changed[i] - expected) > 1e-8 * (1. + std::abs(expected)))
          throw std::runtime_error("Relaxation length changed the wrong equation");
      }
    }
  }
  for (double bad : {0., -.1}) {
    parameters->set_value(key, bad);
    bool rejected = false;
    try {
      model.evaluate();
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    if (!rejected) throw std::runtime_error("Invalid relaxation length accepted");
  }
}
