#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ocd_steering_actuator_pt1_cpp/steering_actuator_pt1_model.hpp"
using Model = tam::ocd::steering_actuator::PT1SteeringActuatorModel;

int main()
{
  Model model;
  auto parameters = model.get_param_manager();
  parameters->set_value("steering_actuator.use_damped_pt2", true);
  parameters->set_value("steering_actuator.natural_frequency_hz", 4.);
  parameters->set_value("steering_actuator.damping_ratio", .35);
  parameters->set_value("steering_actuator.angle_max_rad", .5);
  parameters->set_value("steering_actuator.angle_rate_max_radps", 10.);
  model.set_driver_input({.2});
  auto state = Model::StateVectorType::Zero().eval();
  double peak = 0.;
  constexpr double dt = .0005;
  for (int i = 0; i < 2000; ++i) {
    model.set_x_vec(state);
    model.evaluate();
    const auto derivative = model.get_x_dot_vec();
    if (std::abs(derivative[0]) > 10. + 1e-12)
      throw std::runtime_error("PT2 exceeded slew limit");
    state += dt * derivative;
    peak = std::max(peak, state[0]);
  }
  if (!(peak > .21 && peak < .5 && std::abs(state[0] - .2) < .001))
    throw std::runtime_error("PT2 overshoot or settling incorrect");
}
