#include <cmath>
#include <stdexcept>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace states = tam::ocd::vehicle_dynamics::double_track::States;

double scalar(const tam::tsl::data_map_t & data, const std::string & key)
{
  return std::visit(
    [](const auto & v) -> double {
      using T = std::decay_t<decltype(v)>;
      if constexpr (std::is_arithmetic_v<T>) return v;
      else if constexpr (std::is_pointer_v<T>) return *v;
      else return NAN;
    }, data.at(key));
}

int main()
{
  Model model;
  const auto parameters = model.get_param_manager();
  const std::string base = "vehicle_dynamics_double_track.suspension.damper_curve.";
  parameters->set_value(base + "enabled.front", 1.);
  parameters->set_value(base + "wheel_damper_motion_ratio.front", 1.);
  parameters->set_value(base + "speed_mps", std::vector<double>{0., .1, .3});
  parameters->set_value(base + "bump_force_N.front", std::vector<double>{0., 150., 400.});
  parameters->set_value(base + "rebound_force_N.front", std::vector<double>{0., 100., 300.});
  auto state = Model::StateVectorType::Zero().eval();
  for (double speed : {-.4, -.2, -.05, .05, .2, .4}) {
    state[states::v_z_mps] = speed;
    model.set_x_vec(state);
    model.evaluate();
    const auto data = model.get_logger()->get_data();
    const double shaft = scalar(data, "imr/suspension_damper_compression_speed_mps/front_left");
    const double force = scalar(data, "imr/suspension_damper_force_N/front_left");
    const double extrapolated = scalar(data, "imr/damper_curve_speed_extrapolated/front_left");
    if (!std::isfinite(force) || force * shaft < 0. ||
        (std::abs(shaft) > .3) != (extrapolated > .5))
      throw std::runtime_error("Damper response is nonpassive or extrapolation unreported");
  }
  parameters->set_value(base + "bump_force_N.front", std::vector<double>{0., -10., 400.});
  bool rejected = false;
  try {
    model.evaluate();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Negative damper force curve accepted");
}
