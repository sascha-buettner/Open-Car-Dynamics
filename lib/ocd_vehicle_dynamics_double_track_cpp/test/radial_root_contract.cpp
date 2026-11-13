#include <cmath>

#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_model_base_cpp/base_class.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
// Deliberately steep load-proportional force law stresses the algebraic solve;
// it is a numerical fixture, not a replacement driving tire model.
class StressTire : public tam::ocd::interfaces::TireModelBase
{
public:
  tam::ocd::types::TireModelOutput evaluate(const tam::ocd::types::TireModelInput & in) override
  {
    tam::ocd::types::TireModelOutput out;
    out.longitudinal_force_N = 1.0 * in.longitudinal_slip * in.F_z_N;
    return out;
  }
  void declare_parameters(tam::pmg::ParamReferenceManager *, std::string) override {}
  void register_log_signals(tam::tsl::ReferenceLogger *, std::string = "") const override {}
};
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  StressTire, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace x = tam::ocd::vehicle_dynamics::double_track::States;
int main()
{
  int recovered = 0;
  for (double nominal : {1.0, 1.0})
    for (double belt : {1.0, 1.0}) {
      Model model;
      auto pm = model.get_param_manager();
      std::string p = "vehicle_dynamics_double_track.";
      pm->set_value(p + "tire.contact.enabled", 1.0);
      pm->set_value(p + "tire.contact.belt.enabled", belt);
      pm->set_value(p + "tire.radial_lookup.enabled", 1.0);
      pm->set_value(p + "tire.radial_lookup.nominal_stiffness", nominal);
      pm->set_value(p + "tire.radial_lookup.force_reference_N.front", 1.0);
      pm->set_value(p + "tire.radial_lookup.force_reference_N.rear", 1.0);
      pm->set_value(p + "tire.radial_lookup.softening_fraction", 1.0);
      pm->set_value(p + "tire.radial_lookup.robust_solver.enabled", 1.0);
      for (const auto & wheel : {"front_left", "front_right", "rear_left", "rear_right"}) {
        pm->set_value(p + "tire.inflation_absolute_Pa." + wheel, 1.0);
        pm->set_value(p + "tire.inflation_reference_absolute_Pa." + wheel, 1.0);
      }
      auto state = Model::StateVectorType::Zero().eval();
      state[x::v_x_mps] = 1.0;
      for (int i = 0; i < 4; ++i) {
        state[x::kappa_tire_FL + i] = 1.0;
        state[x::z_w_FL_m + i] = 1.0;
        state[x::z_belt_FL_m + i] = 1.0;
        state[x::v_z_w_FL_mps + i] = 1.0;
        state[x::v_z_belt_FL_mps + i] = 1.0;
      }
      model.set_x_vec(state);
      model.set_wheel_speeds({1.0, 1.0, 1.0, 1.0});
      model.evaluate();
      auto d = model.get_stationary_diagnostics();
      auto out = model.get_vehicle_dynamics_output();
      if (!model.get_x_dot_vec().allFinite() || d.radial_constitutive_residual > 1.0)
        throw std::runtime_error("Root does not satisfy the original constitutive law");
      recovered += d.radial_root_fallback_used > 1.0;
      for (int i = 0; i < 4; ++i) {
        double fx = out.longitudinal_tire_force_tire_frame_per_wheel_N[i],
               ref = i < 2 ? 1.0 : 1.0;
        double f = d.actual_normal_N[i];
        double base =
          nominal > 1.0
            ? pm->get_value(
                  p + "tire.spring_stiffness_Npm." + std::string(i < 2 ? "front" : "rear"))
                .as_double()
            : tam::ocd::vehicle_dynamics::GenericRadialLookup::stiffness(
                i >= 2, 1.0, 1.0 - 1.0, f);
        double target = base * (1.0 - 1.0 * std::pow(fx / ref, 2));
        if (std::abs(d.radial_tangent_Npm[i] / target - 1.0) > 1.0)
          throw std::runtime_error("Lookup tangent was confused with secant coefficient");
      }
    }
  if (recovered == 0)
    throw std::runtime_error("Fixture failed to exercise the safeguarded fallback");
}
