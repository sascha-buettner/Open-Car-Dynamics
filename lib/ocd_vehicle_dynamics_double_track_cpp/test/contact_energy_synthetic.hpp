#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
#include "synthetic_vehicle.hpp"

namespace contact_energy_synthetic
{
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
namespace st = tam::ocd::vehicle_dynamics::double_track::States;

struct Surface : tam::ocd::types::RoadSurface
{
  tam::ocd::types::RoadSurfaceSample sample(double x, double y) const override
  {
    return {.005 * std::sin(2. * x) + .004 * std::cos(3. * y),
            .010 * std::cos(2. * x), -.012 * std::sin(3. * y), true};
  }
};

double value(const tam::tsl::data_map_t & data, const std::string & key)
{
  return std::visit(
    [](const auto & v) -> double {
      using T = std::decay_t<decltype(v)>;
      if constexpr (std::is_arithmetic_v<T>) return v;
      else if constexpr (std::is_pointer_v<T>) return *v;
      else throw std::runtime_error("Non-scalar contact log");
    }, data.at(key));
}

void run(bool five_point, bool envelope)
{
  Model model;
  if (five_point) {
    synthetic_vehicle::contact(model);
    model.get_param_manager()->set_value(
      "vehicle_dynamics_double_track.tire.contact.enveloping_radius_m.front",
      envelope ? .05 : 0.);
    model.get_param_manager()->set_value(
      "vehicle_dynamics_double_track.tire.contact.enveloping_radius_m.rear",
      envelope ? .05 : 0.);
  }
  model.set_road_surface(std::make_shared<Surface>());
  auto state = Model::StateVectorType::Zero().eval();
  state[st::x_m] = .3;
  state[st::y_m] = .2;
  state[st::psi_rad] = .15;
  const auto evaluate = [&](const auto & s) {
    model.set_x_vec(s);
    model.evaluate();
    return model.get_logger()->get_data();
  };
  const auto data = evaluate(state);
  double energy = 0., force_x = 0., force_y = 0.;
  for (const char * wheel : {"front_left", "front_right", "rear_left", "rear_right"}) {
    const std::string name(wheel);
    energy += value(data, "imr/contact_elastic_energy_J/" + name);
    force_x += value(data, "imr/contact_force_x_N/" + name);
    force_y += value(data, "imr/contact_force_y_N/" + name);
  }
  if (!(std::isfinite(energy) && energy > 0.))
    throw std::runtime_error("Contact elastic energy is missing");
  const double c = std::cos(state[st::psi_rad]), s = std::sin(state[st::psi_rad]);
  const double global_force[2] = {c * force_x - s * force_y, s * force_x + c * force_y};
  for (int axis = 0; axis < 2; ++axis) {
    const int index = axis == 0 ? st::x_m : st::y_m;
    const double step = 1e-5;
    auto plus = state, minus = state;
    plus[index] += step;
    minus[index] -= step;
    const auto total_energy = [&](const auto & q) {
      const auto d = evaluate(q);
      double sum = 0.;
      for (const char * wheel : {"front_left", "front_right", "rear_left", "rear_right"})
        sum += value(d, std::string("imr/contact_elastic_energy_J/") + wheel);
      return sum;
    };
    const double gradient_force = -(total_energy(plus) - total_energy(minus)) / (2. * step);
    const double tolerance = .02 * std::max(1., std::abs(gradient_force));
    if (!std::isfinite(gradient_force) ||
        std::abs(gradient_force - global_force[axis]) > tolerance)
      throw std::runtime_error("Contact force disagrees with elastic energy gradient");
  }
}
}  // namespace contact_energy_synthetic
