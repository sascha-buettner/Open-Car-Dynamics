#include <cmath>
#include <stdexcept>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
#include "synthetic_vehicle.hpp"
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
  synthetic_vehicle::kinematics(model);
  auto state = Model::StateVectorType::Zero().eval();
  state[states::v_x_mps] = 20.;
  state[states::z_w_FL_m] = .02;
  state[states::v_z_w_FL_mps] = .1;
  model.set_x_vec(state);
  const auto parameters = model.get_param_manager();
  const std::string key =
    "vehicle_dynamics_double_track.suspension.kinematics.material_velocity_enabled";
  parameters->set_value(key, 0.);
  model.evaluate();
  const auto legacy = model.get_logger()->get_data();
  const double old_fl = scalar(legacy,
    "imr/velocity_wheel_over_ground_mps/front_left/longitudinal");
  const double old_fr = scalar(legacy,
    "imr/velocity_wheel_over_ground_mps/front_right/longitudinal");
  parameters->set_value(key, 1.);
  model.evaluate();
  const auto corrected = model.get_logger()->get_data();
  const double new_fl = scalar(corrected,
    "imr/velocity_wheel_over_ground_mps/front_left/longitudinal");
  const double new_fr = scalar(corrected,
    "imr/velocity_wheel_over_ground_mps/front_right/longitudinal");
  if (!model.get_x_dot_vec().allFinite() ||
      std::abs((new_fl - old_fl) - .01) > 1e-9 ||
      std::abs(new_fr - old_fr) > 1e-9)
    throw std::runtime_error("Material Jacobian velocity does not follow local wheel travel");
}
