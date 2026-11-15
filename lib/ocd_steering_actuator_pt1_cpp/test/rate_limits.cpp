#include <cmath>
#include <stdexcept>

#include "ocd_steering_actuator_pt1_cpp/steering_actuator_pt1_model.hpp"
using Model = tam::ocd::steering_actuator::PT1SteeringActuatorModel;

int main()
{
  Model model;
  auto parameters = model.get_param_manager();
  parameters->set_value("steering_actuator.T_PT1", .1);
  parameters->set_value("steering_actuator.angle_max_rad", .3);
  parameters->set_value("steering_actuator.angle_rate_max_radps", .5);
  const auto derivative = [&](double angle, double command) {
    auto state = Model::StateVectorType::Zero().eval();
    state[0] = angle;
    model.set_x_vec(state);
    model.set_driver_input({command});
    model.evaluate();
    return model.get_x_dot_vec()[0];
  };
  if (std::abs(derivative(0., .2) - .5) > 1e-12 ||
      std::abs(derivative(0., -.2) + .5) > 1e-12)
    throw std::runtime_error("Legacy symmetric slew limit changed");
  parameters->set_value("steering_actuator.angle_rate_rise_max_radps", .25);
  parameters->set_value("steering_actuator.angle_rate_fall_max_radps", .1);
  if (std::abs(derivative(0., .2) - .25) > 1e-12 ||
      std::abs(derivative(0., -.2) + .1) > 1e-12)
    throw std::runtime_error("Directional slew limit incorrect");
  if (derivative(.3, .4) != 0. || derivative(-.3, -.4) != 0.)
    throw std::runtime_error("Steering angle stop permits outward motion");
  parameters->set_value("steering_actuator.T_PT1", 0.);
  bool rejected = false;
  try {
    derivative(0., .1);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Invalid steering pole accepted");
}
