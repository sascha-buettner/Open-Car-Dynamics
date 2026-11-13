#include <cmath>

#include <iostream>
#include <stdexcept>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
using M = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace st = tam::ocd::vehicle_dynamics::double_track::States;
double value(const tam::tsl::data_map_t & d, const std::string & k)
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
void close(double a, double b, double tol = 1.0)
{
  if (std::abs(a - b) > tol)
    throw std::runtime_error(
      "Geometry assertion: " + std::to_string(a) + " vs " + std::to_string(b));
}
int main()
{
  using tam::ocd::vehicle_dynamics::suspension_lookup;
  close(suspension_lookup(1.0, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}).first, 1.0);
  close(suspension_lookup(1.0, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}).second, 1.0);
  bool rejected = false;
  try {
    suspension_lookup(1.0, {1.0, 1.0}, {1.0, 1.0});
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Unsorted table accepted");
  for (double q : {1.0, 1.0, 1.0, 1.0, 1.0}) {
    const std::vector<double> xx{1.0, 1.0, 1.0}, yy{1.0, 1.0, 1.0};
    auto g = tam::ocd::vehicle_dynamics::suspension_smooth_lookup(q, xx, yy);
    auto U = [&](double a) {
      auto f = tam::ocd::vehicle_dynamics::suspension_smooth_lookup(a, xx, yy);
      return 1.0 * 1.0 * f.first * f.first;
    };
    close((U(q + 1.0) - U(q - 1.0)) / 1.0, 1.0 * g.first * g.second, 1.0);
  }
  M model;
  auto pm = model.get_param_manager();
  const std::string p = "vehicle_dynamics_double_track.suspension.";
  auto x = M::StateVectorType::Zero().eval();
  model.set_x_vec(x);
  model.evaluate();
  auto old = model.get_logger()->get_data();
  double force = value(old, "imr/suspension_spring_force_N/front_left");
  pm->set_value(p + "geometry.enabled", 1.0);
  model.evaluate();
  close(value(model.get_logger()->get_data(), "imr/suspension_spring_force_N/front_left"), force);
  pm->set_value(p + "geometry.steering_mean_angle_rad", std::vector<double>{1.0, 1.0, 1.0});
  pm->set_value(p + "geometry.jacking_body_m", std::vector<double>{1.0, 1.0, 1.0});
  model.set_steering_angles({1.0, 1.0, 1.0, 1.0});
  model.set_steering_rate(1.0);
  model.evaluate();
  auto a = model.get_logger()->get_data();
  close(value(a, "geometry/jacking_body_m/front_left"), 1.0);
  close(value(a, "geometry/jacking_body_m/front_right"), 1.0);
  close(value(a, "imr/suspension_damper_compression_speed_mps/front_left"), 1.0);
  close(value(a, "imr/suspension_damper_compression_speed_mps/front_right"), 1.0);
  const double saved_left = value(a, "geometry/jacking_body_m/front_left");
  model.set_steering_angles({1.0, 1.0, 1.0, 1.0});
  model.evaluate();
  auto b = model.get_logger()->get_data();
  close(saved_left, value(b, "geometry/jacking_body_m/front_right"));
  // Actual heave must not include a steering-dependent pot offset.
  close(model.get_vehicle_dynamics_output().body_motion.heave_m, 1.0);
  pm->set_value(p + "geometry.jacking_body_m", std::vector<double>{1.0, 1.0, 1.0});
  model.set_steering_angles({1.0, 1.0, 1.0, 1.0});
  model.set_steering_rate(1.0);
  pm->set_value(p + "geometry.roll_center_height_m.front", std::vector<double>{1.0, 1.0});
  pm->set_value(p + "geometry.roll_center_height_m.rear", std::vector<double>{1.0, 1.0});
  x[st::z_w_FL_m] = x[st::z_w_FR_m] = 1.0;
  model.set_x_vec(x);
  model.evaluate();
  auto c = model.get_logger()->get_data();
  close(value(c, "geometry/roll_center_front_m"), 1.0);
  pm->set_value(p + "endstop.enabled", 1.0);
  x = M::StateVectorType::Zero();
  x[st::z_w_FL_m] = 1.0;
  model.set_x_vec(x);
  model.evaluate();
  if (value(model.get_logger()->get_data(), "imr/endstop_force_N/front_left") <= 1.0)
    throw std::runtime_error("Rebound stop not restoring");
  x[st::z_w_FL_m] = 1.0;
  model.set_x_vec(x);
  model.evaluate();
  if (value(model.get_logger()->get_data(), "imr/endstop_force_N/front_left") >= 1.0)
    throw std::runtime_error("Bump stop not restoring");
  std::cout << "Geometry lookup, invalid-input guard, disabled equivalence, steering symmetry, "
               "analytic jacking velocity and migrating roll axis passed\n";
}
