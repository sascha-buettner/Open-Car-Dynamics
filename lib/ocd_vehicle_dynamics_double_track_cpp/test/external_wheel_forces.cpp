#include <cmath>

#include <iostream>
#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace st = tam::ocd::vehicle_dynamics::double_track::States;
int main()
{
  for (bool belt : {false, true}) {
    Model m;
    auto p = m.get_param_manager();
    p->set_value("vehicle_dynamics_double_track.mass_wheel_kg.front", 1.0);
    p->set_value("vehicle_dynamics_double_track.tire.contact.enabled", 1.0);
    p->set_value("vehicle_dynamics_double_track.tire.contact.belt.enabled", double(belt));
    p->set_value("vehicle_dynamics_double_track.tire.contact.belt.mass_kg.front", 1.0);
    m.set_x_vec(Model::StateVectorType::Zero().eval());
    m.evaluate();
    auto rest = m.get_x_dot_vec();
    {
      const auto before = m.get_tire_contact_diagnostics();
      auto changed = m.get_x_vec();
      changed[belt ? st::z_belt_FL_m : st::z_w_FL_m] = 1.0;
      m.set_x_vec(changed);
      m.evaluate();
      const auto after = m.get_tire_contact_diagnostics();
      if (!(after.tire_spring_force_N.front_left > before.tire_spring_force_N.front_left))
        throw std::runtime_error(
          "Live tire diagnostics failed to update after upward wheel/belt motion");
      if (belt && std::abs(after.belt_displacement_m.front_left - 1.0) > 1.0)
        throw std::runtime_error("Belt diagnostics returned a stale position");
      m.set_x_vec(Model::StateVectorType::Zero().eval());
      m.evaluate();
    }
    for (bool tire : {false, true}) {
      tam::ocd::types::ExternalInfluences e{};
      if (tire)
        e.tire_force_z_N.front_left = 1.0;
      else
        e.rim_force_z_N.front_left = 1.0;
      m.set_external_influences(e);
      m.evaluate();
      auto delta = (m.get_x_dot_vec() - rest).eval();
      double momentum_rate =
        (belt ? 1.0 : 1.0) * delta[st::v_z_w_FL_mps] + (belt ? 1.0 : 1.0) * delta[st::v_z_belt_FL_mps];
      if (std::abs(momentum_rate - 1.0) > 1.0)
        throw std::runtime_error("External force violates unsprung momentum balance");
      if (
        std::abs(delta[st::v_z_mps]) + std::abs(delta[st::phi_dot_radps]) +
          std::abs(delta[st::theta_dot_radps]) >
        1.0)
        throw std::runtime_error("Unsprung force also injected as an immediate body wrench");
      if (
        std::abs(delta[st::v_z_w_FR_mps]) + std::abs(delta[st::v_z_w_RL_mps]) +
          std::abs(delta[st::v_z_w_RR_mps]) >
        1.0)
        throw std::runtime_error("Force applied to wrong corner");
      if (belt && std::abs(delta[tire ? st::v_z_w_FL_mps : st::v_z_belt_FL_mps]) > 1.0)
        throw std::runtime_error("Tire/rim force ports collapsed despite separate masses");
    }
  }
  std::cout
    << "Tire/rim force location, unsprung momentum and zero duplicate body wrench verified\n";
}
