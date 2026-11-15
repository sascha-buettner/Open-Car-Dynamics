#include <cmath>
#include <stdexcept>

#include "ocd_steering_actuator_pt1_cpp/steering_actuator_pt1_model.hpp"
using Model = tam::ocd::steering_actuator::PT1SteeringActuatorModel;

int main()
{
  Model model;
  auto parameters = model.get_param_manager();
  auto state = Model::StateVectorType::Zero().eval();
  state[0] = .1;
  model.set_x_vec(state);
  model.evaluate();
  auto plain = model.get_steering_angles();
  if (plain.front_left != .1 || plain.front_right != .1 ||
      plain.rear_left != 0. || plain.rear_right != 0.)
    throw std::runtime_error("Unconfigured steering is not parallel");

  parameters->set_value("steering_actuator.ackermann.track_m", 1.5);
  parameters->set_value("steering_actuator.ackermann.wheelbase_m", 3.);
  parameters->set_value("steering_actuator.ackermann.mean_angle_rad",
                        std::vector<double>{0., .1, .3});
  parameters->set_value("steering_actuator.ackermann.fraction",
                        std::vector<double>{0., .5, 1.});
  for (double angle : {-.2, -.1, .1, .2}) {
    state[0] = angle;
    model.set_x_vec(state);
    model.evaluate();
    const auto wheel = model.get_steering_angles();
    if (std::abs((wheel.front_left + wheel.front_right) / 2. - angle) > 1e-12)
      throw std::runtime_error("Ackermann changed mean steer angle");
    const double inner = std::max(std::abs(wheel.front_left), std::abs(wheel.front_right));
    const double outer = std::min(std::abs(wheel.front_left), std::abs(wheel.front_right));
    const double actual =
      (1. / std::tan(outer) - 1. / std::tan(inner)) / (1.5 / 3.);
    const double expected = .5 + (std::abs(angle) - .1) * 2.5;
    if (std::abs(actual - expected) > 1e-6)
      throw std::runtime_error("Ackermann fraction incorrect");
    if ((angle > 0.) != (std::abs(wheel.front_left) > std::abs(wheel.front_right)))
      throw std::runtime_error("Ackermann inner wheel mirrored incorrectly");
  }
  parameters->set_value("steering_actuator.ackermann.fraction",
                        std::vector<double>{0., .5, 1.2});
  bool rejected = false;
  try {
    model.get_steering_angles();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Invalid Ackermann fraction accepted");
}
