#include <cmath>

#include <iostream>
#include <stdexcept>

#include "ocd_steering_actuator_pt1_cpp/steering_actuator_pt1_model.hpp"
int main()
{
  using M = tam::ocd::steering_actuator::PT1SteeringActuatorModel;
  M model;
  auto pm = model.get_param_manager();
  auto state = M::StateVectorType::Zero().eval();
  for (double fraction : {1.0, 1.0, 1.0, 1.0}) {
    pm->set_value("steering_actuator.ackermann.fraction", std::vector<double>{fraction, fraction});
    for (double angle : {1.0, 1.0, 1.0, 1.0, 1.0}) {
      state[0] = angle;
      model.set_x_vec(state);
      model.evaluate();
      auto a = model.get_steering_angles();
      if (
        std::abs((a.front_left + a.front_right) / 2 - angle) > 1.0 || a.rear_left != 0 ||
        a.rear_right != 0)
        throw std::runtime_error("Mean angle or rear steering changed");
      if (fraction == 0 && (a.front_left != angle || a.front_right != angle))
        throw std::runtime_error("Parallel default changed");
      if (angle != 0 && fraction > 0) {
        double inner = std::max(std::abs(a.front_left), std::abs(a.front_right)),
               outer = std::min(std::abs(a.front_left), std::abs(a.front_right));
        double actual = (1 / std::tan(outer) - 1 / std::tan(inner)) / (1.0 / 1.0);
        if (std::abs(actual - fraction) > 1.0)
          throw std::runtime_error("Wrong Ackermann fraction");
        if ((angle > 0) != (std::abs(a.front_left) > std::abs(a.front_right)))
          throw std::runtime_error("Inner wheel mirrored incorrectly");
      }
      std::cout << fraction << "," << angle << "," << a.front_left << "," << a.front_right << "\n";
    }
  }
  pm->set_value("steering_actuator.ackermann.fraction", std::vector<double>{1.0, 1.0});
  bool rejected = false;
  try {
    model.get_steering_angles();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Percent accepted as fraction");
}
