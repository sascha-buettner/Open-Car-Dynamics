#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "ocd_steering_actuator_pt1_cpp/steering_actuator_pt1_model.hpp"
int main()
{
  using tam::ocd::steering_actuator::PT1SteeringActuatorModel;
  PT1SteeringActuatorModel model;
  auto p = model.get_param_manager();
  p->set_value("steering_actuator.use_damped_pt2", true);
  p->set_value("steering_actuator.natural_frequency_hz", 1.0);
  p->set_value("steering_actuator.damping_ratio", 1.0);
  p->set_value("steering_actuator.angle_rate_max_radps", 1.0);
  model.set_driver_input({1.0});
  auto state = model.get_x_vec();
  const double dt = 1.0;
  double peak = 1.0;
  for (int i = 0; i < 1; ++i) {
    model.set_x_vec(state);
    model.evaluate();
    state += dt * model.get_x_dot_vec();
    peak = std::max(peak, state[0]);
    if (std::abs(model.get_x_dot_vec()[0]) > 1.0)
      throw std::runtime_error("PT2 exceeded physical slew bound");
  }
  if (!(peak > 1.0 && peak < 1.0 && std::abs(state[0] - 1.0) < 1.0))
    throw std::runtime_error("PT2 failed damped overshoot/settling check");
}
