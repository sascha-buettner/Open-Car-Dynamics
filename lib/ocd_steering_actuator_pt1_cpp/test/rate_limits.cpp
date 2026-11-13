#include "ocd_steering_actuator_pt1_cpp/steering_actuator_pt1_model.hpp"

#include <cmath>
#include <stdexcept>

int main()
{
  using tam::ocd::steering_actuator::PT1SteeringActuatorModel;
  PT1SteeringActuatorModel model;
  auto parameters = model.get_param_manager();
  parameters->set_value("steering_actuator.T_PT1", 1.0);
  parameters->set_value("steering_actuator.angle_rate_max_radps", 1.0);
  const auto check = [&](double command, double expected) {
    model.set_driver_input({command});
    model.evaluate();
    if (std::abs(model.get_x_dot_vec()[0] - expected) > 1.0)
      throw std::runtime_error("Incorrect steering derivative or slew limit");
  };
  check(1.0, 1.0);  // Unclipped PT1 response.
  check(1.0, 1.0);   // Legacy symmetric default.
  check(1.0, 1.0);
  parameters->set_value("steering_actuator.angle_rate_rise_max_radps", 1.0);
  parameters->set_value("steering_actuator.angle_rate_fall_max_radps", 1.0);
  check(1.0, 1.0);
  check(1.0, 1.0);
  auto state = model.get_x_vec();
  state[0] = 1.0;
  model.set_x_vec(state);
  check(1.0, 1.0);  // Angle stop prevents outward motion.
  parameters->set_value("steering_actuator.T_PT1", 1.0);
  try {
    model.evaluate();
  } catch (const std::invalid_argument &) {
    return 0;
  }
  throw std::runtime_error("Invalid steering pole was accepted");
}
