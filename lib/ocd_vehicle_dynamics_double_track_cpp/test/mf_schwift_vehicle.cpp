#include <cmath>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/road_contact_frame.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
#include "synthetic_vehicle.hpp"

using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace st = tam::ocd::vehicle_dynamics::double_track::States;
struct Plane : tam::ocd::types::RoadSurface
{
  tam::ocd::types::RoadSurfaceSample sample(double x, double) const override
  {
    return {.05 * x, .05, 0., true};
  }
};

double scalar(const tam::tsl::data_map_t & data, const std::string & key)
{
  return std::visit(
    [](const auto & v) -> double {
      using T = std::decay_t<decltype(v)>;
      if constexpr (std::is_arithmetic_v<T>) return v;
      else if constexpr (std::is_pointer_v<T>) return *v;
      else throw std::runtime_error("Non-scalar ring log");
    }, data.at(key));
}

int main()
{
  Model model;
  synthetic_vehicle::contact(model, true);
  auto pm = model.get_param_manager();
  const std::string base = "vehicle_dynamics_double_track.tire.mfschwift.";
  for (const auto & item : {
         std::pair<const char *, double>{"inertia_transverse_kgm2", .04},
         {"inertia_spin_kgm2", .08}, {"frequency_lateral_hz", 40.},
         {"frequency_yaw_camber_hz", 35.}, {"frequency_windup_hz", 45.},
         {"damping_lateral", .15}, {"damping_yaw_camber", .15},
         {"damping_windup", .15}, {"longitudinal_stiffness_Npm", 100000.}})
    pm->set_value(base + item.first, item.second);
  pm->set_value(base + "enabled", 1.);
  auto state = Model::StateVectorType::Zero().eval();
  model.set_x_vec(state);
  model.evaluate();
  if (!model.get_x_dot_vec().allFinite())
    throw std::runtime_error("Synthetic ring equilibrium is nonfinite");
  state[st::ring_x_FL_m] = .001;
  model.set_x_vec(state);
  model.evaluate();
  const double spring_accel = model.get_x_dot_vec()[st::ring_vx_FL_mps];
  const double energy = scalar(model.get_logger()->get_data(), "mfschwift/ring_energy_J/front_left");
  if (!(spring_accel < 0. && energy > 0.))
    throw std::runtime_error("Ring sidewall does not restore displacement");
  state[st::ring_vx_FL_mps] = .2;
  model.set_x_vec(state);
  model.evaluate();
  const double damping = scalar(model.get_logger()->get_data(), "mfschwift/sidewall_dissipation_W/front_left");
  if (!(std::isfinite(damping) && damping > 0.))
    throw std::runtime_error("Ring sidewall is not dissipative");
  model.set_road_surface(std::make_shared<Plane>());
  state = Model::StateVectorType::Zero().eval();
  state[st::ring_vx_FL_mps] = 8.;
  state[st::ring_omega_y_FL_radps] = 25.;
  state[st::v_z_belt_FL_mps] = 2.;
  state[st::kappa_tire_FL] = .0001;
  model.set_x_vec(state);
  model.evaluate();
  const auto data = model.get_logger()->get_data();
  const double radius = scalar(data, "imr/dynamic_tire_radius/front_left");
  const tam::ocd::vehicle_dynamics::RoadContactFrame frame(.05, 0.);
  const double contact_speed = frame.velocity({8., 0., 2.}).x();
  const double ring_spin_speed = radius * 25. * frame.spin_lever_scale;
  const double expected_slip = std::clamp(
    (ring_spin_speed - contact_speed) / std::max(3., contact_speed), -1., 1.);
  if (std::abs(scalar(data, "imr/tire_longitudinal_slip/front_left") - expected_slip) > 1e-8)
    throw std::runtime_error("Ring slip omits local road-plane velocity");
  const double fx = scalar(data, "imr/tire_forces_tire_frame_N/front_left/longitudinal");
  const double rolling = scalar(data, "imr/rolling_resistance_N/front_left");
  const double sidewall_torque = -scalar(data, "mfschwift/spindle_spin_load_Nm/front_left");
  const double torque = .08 * model.get_x_dot_vec()[st::ring_omega_y_FL_radps] + sidewall_torque;
  const double expected_torque = -fx * radius * frame.spin_lever_scale -
    rolling * radius * std::tanh(25. * radius / .1);
  if (!(std::abs(fx) > 1. && std::abs(torque - expected_torque) < 1e-8))
    throw std::runtime_error("Ring spin torque disagrees with contact-plane rolling lever");
  pm->set_value(base + "inertia_spin_kgm2", 0.);
  bool rejected = false;
  try { model.evaluate(); } catch (const std::invalid_argument &) { rejected = true; }
  if (!rejected) throw std::runtime_error("Zero ring inertia accepted");
}
