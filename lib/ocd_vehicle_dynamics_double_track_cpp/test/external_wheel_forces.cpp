#include <cmath>
#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
#include "synthetic_vehicle.hpp"
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace states = tam::ocd::vehicle_dynamics::double_track::States;

int main()
{
  constexpr double force = 100.;
  for (bool belt : {false, true}) {
    Model model;
    synthetic_vehicle::contact(model, belt);
    model.set_x_vec(Model::StateVectorType::Zero().eval());
    model.evaluate();
    const auto rest = model.get_x_dot_vec().eval();
    for (bool tire_port : {false, true}) {
      tam::ocd::types::ExternalInfluences input{};
      if (tire_port) input.tire_force_z_N.front_left = force;
      else input.rim_force_z_N.front_left = force;
      model.set_external_influences(input);
      model.evaluate();
      const auto change = (model.get_x_dot_vec() - rest).eval();
      const double rim_mass = belt ? 12. : 15.;
      const double belt_mass = belt ? 3. : 0.;
      const double expected_rim = (!belt || !tire_port) ? force / rim_mass : 0.;
      const double expected_belt = belt && tire_port ? force / belt_mass : 0.;
      if (std::abs(change[states::v_z_w_FL_mps] - expected_rim) > 1e-8 ||
          std::abs(change[states::v_z_belt_FL_mps] - expected_belt) > 1e-8 ||
          std::abs(change[states::v_z_mps]) > 1e-8 ||
          std::abs(change[states::v_z_w_FR_mps]) > 1e-8)
        throw std::runtime_error("External force used the wrong tire/rim port");
    }
  }
}
