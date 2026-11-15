#include <cmath>
#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/mf_52.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"

int main()
{
  using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
    tam::ocd::tire_models::MF52, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
  Model model;
  const auto parameters = model.get_param_manager();
  const auto value = [&](const char * name) {
    return parameters->get_value(std::string("vehicle_dynamics_double_track.") + name).as_double();
  };
  // These generic values already existed on public main before the optional mechanisms.
  if (value("mass_vehicle_kg") != 800. || value("wheelbase_m") != 3.2 ||
      value("tire.spring_stiffness_Npm.front") != 250000. ||
      value("z.tire_model.front_left.PCX1") != 1.6 ||
      value("z.tire_model.front_left.FITTYP") != 52.)
    throw std::runtime_error("Public baseline parameters changed");
  for (const char * feature : {
         "suspension.geometry.enabled", "suspension.kinematics.enabled",
         "suspension.endstop.enabled", "suspension.damper_curve.enabled.front",
         "suspension.corner_weights.enabled", "tire.contact.enabled",
         "tire.contact.belt.enabled", "tire.radial_lookup.enabled",
         "tire.mf62_vertical.enabled", "tire.load_memory.fraction",
         "tire.relaxation.longitudinal.recovery_length_m.front"})
    if (value(feature) != 0.) throw std::runtime_error("Optional feature enabled by default");
  auto state = Model::StateVectorType::Zero().eval();
  model.set_x_vec(state);
  model.evaluate();
  if (!model.get_x_dot_vec().allFinite())
    throw std::runtime_error("Public baseline produces nonfinite derivatives");
  state[tam::ocd::vehicle_dynamics::double_track::States::v_x_mps] = 20.;
  model.set_x_vec(state);
  model.set_wheel_speeds({20. / 0.3, 20. / 0.3, 20. / 0.3, 20. / 0.3});
  model.evaluate();
  if (!model.get_x_dot_vec().allFinite())
    throw std::runtime_error("Public baseline fails at rolling speed");
}
