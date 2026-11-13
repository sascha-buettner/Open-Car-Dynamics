// Copyright 2026 Simon Sagmeister

#include "ocd_steering_actuator_pt1_cpp/steering_actuator_pt1_model.hpp"

#include <cmath>

#include <stdexcept>
namespace tam::ocd::steering_actuator
{
PT1SteeringActuatorModel::PT1SteeringActuatorModel()
{
  declare_parameters();
  register_log_signals();
}
void PT1SteeringActuatorModel::evaluate()
{
  double max_angle = p_.angle_max_rad;
  const auto rate = [this](double directional) {
    return directional == -1.0 ? p_.angle_rate_max_radps : directional;
  };
  const double rise = rate(p_.angle_rate_rise_max_radps);
  const double fall = rate(p_.angle_rate_fall_max_radps);
  if (!(std::isfinite(p_.T_PT1) && p_.T_PT1 > 0.0 && std::isfinite(max_angle) && max_angle > 0.0 &&
        std::isfinite(rise) && rise > 0.0 && std::isfinite(fall) && fall > 0.0))
    throw std::invalid_argument("Steering pole, angle and rate limits must be finite and positive");
  double x_dot;
  x_dot_vec_[x::velocity_radps] = 0.0;
  if (p_.use_damped_pt2) {
    if (!(std::isfinite(p_.natural_frequency_hz) && p_.natural_frequency_hz > 0.0 &&
          std::isfinite(p_.damping_ratio) && p_.damping_ratio > 0.0))
      throw std::invalid_argument("Damped steering PT2 frequency and damping must be positive");
    const double omega = 2.0 * std::acos(-1.0) * p_.natural_frequency_hz;
    const double velocity = std::clamp(x_vec_[x::velocity_radps], -fall, rise);
    double acceleration =
      omega * omega * (requested_steering_angle_rad_ - x_vec_[x::position_rad]) -
      2.0 * p_.damping_ratio * omega * velocity;
    // Prevent hidden velocity windup while the physical slew/angle cap is active.
    if (
      (velocity >= rise && acceleration > 0.0) || (velocity <= -fall && acceleration < 0.0) ||
      (x_vec_[x::position_rad] >= max_angle && acceleration > 0.0) ||
      (x_vec_[x::position_rad] <= -max_angle && acceleration < 0.0))
      acceleration = 0.0;
    x_dot_vec_[x::velocity_radps] = acceleration;
    x_dot = velocity;
  } else {
    x_dot = (requested_steering_angle_rad_ - x_vec_[x::position_rad]) / p_.T_PT1;
  }
  // Enforce the maximum rate
  x_dot = std::clamp(x_dot, -fall, rise);
  // Enforce the maximum angle
  if (x_vec_[x::position_rad] <= -max_angle) x_dot = std::max(0.0, x_dot);
  if (x_vec_[x::position_rad] >= max_angle) x_dot = std::min(0.0, x_dot);
  // Set the x_dot vec
  x_dot_vec_[x::position_rad] = x_dot;
  measured_steering_angle_rad_ = x_vec_[x::position_rad];
  actual_steering_angle_rad_ = measured_steering_angle_rad_ + p_.static_offset_rad;
}
void PT1SteeringActuatorModel::set_x_vec(const StateVectorType & x_vec) { x_vec_ = x_vec; }
void PT1SteeringActuatorModel::set_driver_input(const DriverInputType & input)
{
  requested_steering_angle_rad_ = input.steering_angle_rad;
}
void PT1SteeringActuatorModel::set_load(
  const DoublePerWheelType & steering_load_torque_per_wheel_Nm)
{
  steering_load_torque_per_wheel_Nm_ = steering_load_torque_per_wheel_Nm;
}
PT1SteeringActuatorModel::StateVectorType PT1SteeringActuatorModel::get_x_vec() const
{
  return x_vec_;
}
PT1SteeringActuatorModel::StateVectorType PT1SteeringActuatorModel::get_x_dot_vec() const
{
  return x_dot_vec_;
}
tam::tsl::LoggerAccessInterface::SharedPtr PT1SteeringActuatorModel::get_logger() const
{
  return logger_;
}
tam::pmg::MgmtInterface::SharedPtr PT1SteeringActuatorModel::get_param_manager() const
{
  return param_manager_;
}
PT1SteeringActuatorModel::DoublePerWheelType PT1SteeringActuatorModel::get_steering_angles() const
{
  DoublePerWheelType out;
  out.front_left = actual_steering_angle_rad_;
  out.front_right = actual_steering_angle_rad_;
  const auto & knots = p_.ackermann_mean_angle_rad;
  const auto & fractions = p_.ackermann_fraction;
  if (
    knots.size() < 2 || knots.size() != fractions.size() || knots.front() != 0. ||
    !(std::isfinite(p_.ackermann_track_m) && p_.ackermann_track_m > 0.) ||
    !(std::isfinite(p_.ackermann_wheelbase_m) && p_.ackermann_wheelbase_m > 0.))
    throw std::invalid_argument("Invalid Ackermann table/dimensions");
  for (std::size_t i = 0; i < knots.size(); ++i)
    if (
      !(std::isfinite(knots[i]) && std::isfinite(fractions[i])) || fractions[i] < 0. ||
      fractions[i] > 1. || (i && knots[i] <= knots[i - 1]))
      throw std::invalid_argument("Ackermann angles must increase; fractions must be in [0,1]");
  const double mean = std::abs(actual_steering_angle_rad_);
  const auto next = std::upper_bound(knots.begin(), knots.end(), mean);
  const auto hi = std::clamp<std::size_t>(next - knots.begin(), 1, knots.size() - 1);
  const double alpha = std::clamp((mean - knots[hi - 1]) / (knots[hi] - knots[hi - 1]), 0., 1.);
  const double fraction = fractions[hi - 1] + alpha * (fractions[hi] - fractions[hi - 1]);
  if (fraction > 0. && mean > 1e-7) {
    // Preserve the actuator's arithmetic mean wheel angle and feedback exactly.
    // Solve the standard Ackermann identity, rather than adding a toe offset.
    const double target = fraction * p_.ackermann_track_m / p_.ackermann_wheelbase_m;
    double lo = 0., hi_angle = mean * (1. - 1e-12);
    for (int i = 0; i < 48; ++i) {
      const double half_delta = (lo + hi_angle) / 2;
      const double value = 1. / std::tan(mean - half_delta) - 1. / std::tan(mean + half_delta);
      if (value < target)
        lo = half_delta;
      else
        hi_angle = half_delta;
    }
    const double half_delta = (lo + hi_angle) / 2;
    const double direction = std::copysign(1., actual_steering_angle_rad_);
    out.front_left = direction * (mean + direction * half_delta);
    out.front_right = direction * (mean - direction * half_delta);
  }
  out.rear_left = 0;
  out.rear_right = 0;
  return out;
}
PT1SteeringActuatorModel::FeedbackType PT1SteeringActuatorModel::get_feedback() const
{
  FeedbackType out;
  out.steering_angle_rad = measured_steering_angle_rad_;
  return out;
}
void PT1SteeringActuatorModel::declare_parameters()
{
  using pt = tam::pmg::ParameterType;
  param_manager_->declare_parameter(
    "steering_actuator.use_damped_pt2", &p_.use_damped_pt2, false, pt::BOOL,
    "Use a damped second-order physical actuator instead of the legacy PT1");
  param_manager_->declare_parameter(
    "steering_actuator.natural_frequency_hz", &p_.natural_frequency_hz, 1.0, pt::DOUBLE,
    "Damped PT2 natural frequency in Hz");
  param_manager_->declare_parameter(
    "steering_actuator.damping_ratio", &p_.damping_ratio, 1.0, pt::DOUBLE,
    "Damped PT2 dimensionless damping ratio; below one permits overshoot");
  param_manager_->declare_parameter(
    "steering_actuator.ackermann.track_m", &p_.ackermann_track_m, 1.0, pt::DOUBLE,
    "Ground contact track [m]");
  param_manager_->declare_parameter(
    "steering_actuator.ackermann.wheelbase_m", &p_.ackermann_wheelbase_m, 1.0, pt::DOUBLE,
    "Wheelbase [m]");
  param_manager_->declare_parameter(
    "steering_actuator.ackermann.mean_angle_rad", &p_.ackermann_mean_angle_rad,
    std::vector<double>{1.0, 1.0}, pt::DOUBLE_ARRAY, "Absolute mean tire-angle knots [rad]");
  param_manager_->declare_parameter(
    "steering_actuator.ackermann.fraction", &p_.ackermann_fraction, std::vector<double>{1.0, 1.0},
    pt::DOUBLE_ARRAY, "Standard Ackermann fraction, not percent");
  param_manager_->declare_parameter("steering_actuator.T_PT1", &p_.T_PT1, 0.06, pt::DOUBLE, "");
  param_manager_->declare_parameter(
    "steering_actuator.angle_max_rad", &p_.angle_max_rad, 0.3, pt::DOUBLE, "");
  param_manager_->declare_parameter(
    "steering_actuator.angle_rate_max_radps", &p_.angle_rate_max_radps, 0.5, pt::DOUBLE, "");
  param_manager_->declare_parameter(
    "steering_actuator.angle_rate_rise_max_radps", &p_.angle_rate_rise_max_radps, 1.0, pt::DOUBLE,
    "Positive tire-angle slew limit; -1 inherits angle_rate_max_radps");
  param_manager_->declare_parameter(
    "steering_actuator.angle_rate_fall_max_radps", &p_.angle_rate_fall_max_radps, 1.0, pt::DOUBLE,
    "Negative tire-angle slew magnitude; -1 inherits angle_rate_max_radps");
  param_manager_->declare_parameter(
    "steering_actuator.static_offset_rad", &p_.static_offset_rad, 0.0, pt::DOUBLE,
    "'Oval steering'");
}
void PT1SteeringActuatorModel::register_log_signals()
{
  logger_->log("target_steering_angle_rad", &requested_steering_angle_rad_);
  logger_->log("measured_steering_angle_rad", &measured_steering_angle_rad_);
  logger_->log("actual_steering_angle_rad", &actual_steering_angle_rad_);
  logger_->log("steering_angle_rate_radps", &x_dot_vec_[x::position_rad]);
}
}  // namespace tam::ocd::steering_actuator
