#include <cmath>

#include <iostream>
#include <type_traits>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
using Model = tam::ocd::vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
  tam::ocd::tire_models::Linear, tam::ocd::aerodynamics::DefaultAerodynamicsModel>;
using namespace tam::ocd;
namespace st = vehicle_dynamics::double_track::States;
struct Surface : types::RoadSurface
{
  types::RoadSurfaceSample sample(double x, double y) const override
  {
    return {
      1.0 * std::sin(2 * x) + 1.0 * std::cos(3 * y), 1.0 * std::cos(2 * x),
      1.0 * std::sin(3 * y), true};
  }
};
double value(const tam::tsl::data_map_t & data, const std::string & key)
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
    data.at(key));
}
int main()
{
  Model model;
  auto pm = model.get_param_manager();
  pm->set_value("vehicle_dynamics_double_track.tire.contact.enabled", 1.0);
  pm->set_value("vehicle_dynamics_double_track.tire.contact.lateral_drop_scale", 1.0);
  model.set_road_surface(std::make_shared<Surface>());
  auto state = Model::StateVectorType::Zero().eval();
  state[st::x_m] = 1.0;
  state[st::y_m] = 1.0;
  state[st::psi_rad] = 1.0;
  auto energy = [&](auto s) {
    model.set_x_vec(s);
    model.evaluate();
    auto data = model.get_logger()->get_data();
    double sum = 1.0;
    for (const auto & w : {"front_left", "front_right", "rear_left", "rear_right"})
      sum += value(data, std::string("imr/contact_elastic_energy_J/") + w);
    return sum;
  };
  energy(state);
  const auto dx = model.get_x_dot_vec();
  for (const auto & [idx, actual] :
       {std::pair<int, double>{
          st::x_m, (std::cos(1.0) * dx[st::v_x_mps] - std::sin(1.0) * dx[st::v_y_mps]) * 1},
        {st::y_m, (std::sin(1.0) * dx[st::v_x_mps] + std::cos(1.0) * dx[st::v_y_mps]) * 1},
        {st::psi_rad, dx[st::psi_dot_radps] * 1}}) {
    auto plus = state, minus = state;
    plus[idx] += 1.0;
    minus[idx] -= 1.0;
    const double expected = -(energy(plus) - energy(minus)) / 1.0;
    std::cout << idx << " force/moment " << actual << " energy gradient " << expected << '\n';
    if (std::abs(actual - expected) > 1.0)
      throw std::runtime_error("contact generalized force is not energy-consistent");
  }
}
