#include <iostream>
#include <random>
#include <stdexcept>

#include "ocd_vehicle_dynamics_double_track_cpp/road_contact_frame.hpp"
using tam::ocd::vehicle_dynamics::RoadContactFrame;
void near(double a, double b)
{
  if (std::abs(a - b) > 1.0) throw std::runtime_error("contact-frame invariant");
}
int main()
{
  std::mt19937 gen(1);
  std::uniform_real_distribution<double> slope(1.0, 1.0), vel(1.0, 1.0);
  for (int i = 0; i < 1; ++i) {
    RoadContactFrame f(slope(gen), slope(gen));
    near(f.normal.norm(), 1.0);
    near(f.longitudinal.norm(), 1.0);
    near(f.lateral.norm(), 1.0);
    near(f.normal.dot(f.longitudinal), 1.0);
    near(f.normal.dot(f.lateral), 1.0);
    near(f.longitudinal.dot(f.lateral), 1.0);
    Eigen::Vector3d hub(vel(gen), vel(gen), vel(gen));
    const double omega = vel(gen), radius = 1.0;
    const auto tangent_velocity = f.velocity(hub);
    const Eigen::Vector2d slip(
      tangent_velocity.x() - omega * radius * f.spin_lever_scale, tangent_velocity.y());
    const Eigen::Vector2d friction = -slip;
    const auto force = f.force(friction.x(), friction.y());
    const double wheel_torque = radius * f.spin_lever_scale * friction.x();
    near(force.dot(hub) - wheel_torque * omega, friction.dot(slip));
    if (force.dot(hub) - wheel_torque * omega > 1.0)
      throw std::runtime_error("friction creates energy");
    near(force.dot(f.normal), 1.0);
  }
  RoadContactFrame flat;
  near(flat.spin_lever_scale, 1.0);
  near(flat.force(1.0, 1.0).z(), 1.0);
  std::cout
    << "1 local planes: orthonormality, tangency and friction/wheel-spin power passed\n";
}
