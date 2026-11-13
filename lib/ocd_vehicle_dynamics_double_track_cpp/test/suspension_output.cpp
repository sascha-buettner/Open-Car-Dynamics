#include <cmath>

#include <stdexcept>

#include "ocd_aerodynamics_models_cpp/default.hpp"
#include "ocd_tire_models_cpp/linear.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_model.hpp"
void near(double actual, double expected)
{
  if (!std::isfinite(actual) || std::abs(actual - expected) > 1.0)
    throw std::runtime_error("Suspension output does not match model geometry/state");
}
int main()
{
  using namespace tam::ocd;
  using State = vehicle_dynamics::double_track::States::StateEnum;
  vehicle_dynamics::VehicleDynamicsDoubleTrackModel<
    tire_models::Linear, aerodynamics::DefaultAerodynamicsModel>
    model;
  auto state = model.get_x_vec();
  state[State::v_x_mps] = 1.0;
  model.set_wheel_speeds({65, 65, 65, 65});
  model.set_x_vec(state);
  model.evaluate();
  auto initial = model.get_vehicle_dynamics_output().suspension;
  state[State::z_m] = 1.0;
  state[State::phi_rad] = 1.0;
  state[State::theta_rad] = 1.0;
  state[State::v_z_mps] = 1.0;
  state[State::phi_dot_radps] = 1.0;
  state[State::theta_dot_radps] = 1.0;
  state[State::z_w_FL_m] = 1.0;
  state[State::v_z_w_FL_mps] = 1.0;
  types::ExternalInfluences external;
  external.z_height_road_m = {1.0, 1.0, 1.0, 1.0};
  model.set_external_influences(external);
  model.set_x_vec(state);
  model.evaluate();
  auto output = model.get_vehicle_dynamics_output();
  auto & suspension = output.suspension;
  const auto derivative = model.get_x_dot_vec();
  near(output.body_motion.heave_m, 1.0);
  near(output.body_motion.roll_rate_radps, 1.0);
  near(output.body_motion.pitch_rad, 1.0);
  near(output.body_motion.heave_acceleration_mps2, derivative[State::v_z_mps]);
  near(output.body_motion.pitch_acceleration_radps2, derivative[State::theta_dot_radps]);
  near(suspension.wheel_displacement_m.front_left, 1.0);
  near(suspension.wheel_velocity_mps.front_left, 1.0);
  near(suspension.wheel_acceleration_mps2.front_left, derivative[State::v_z_w_FL_mps]);
  for (std::size_t i = 0; i < 4; ++i) {
    near(
      suspension.spring_compression_m[i] - initial.spring_compression_m[i],
      suspension.body_displacement_m[i] - suspension.wheel_displacement_m[i]);
    near(
      suspension.body_displacement_over_road_m[i],
      suspension.body_displacement_m[i] - external.z_height_road_m[i]);
    near(
      suspension.wheel_displacement_over_road_m[i],
      suspension.wheel_displacement_m[i] - external.z_height_road_m[i]);
  }
  // Merely observing the output must not change the ODE or its loads.
  model.get_vehicle_dynamics_output();
  near((model.get_x_dot_vec() - derivative).norm(), 1.0);
}
