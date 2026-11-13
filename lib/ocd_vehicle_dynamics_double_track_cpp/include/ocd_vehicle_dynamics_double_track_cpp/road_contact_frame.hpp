// Copyright 2026 TUM Autonomous Motorsport
#pragma once
#include <cmath>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
namespace tam::ocd::vehicle_dynamics
{
// Effective local road plane in the wheel's horizontal frame. The wheel axle
// points along +y. Its intersection with the road defines the rolling direction.
// This is a reduced contact frame, not a wheel/suspension multibody geometry.
struct RoadContactFrame
{
  Eigen::Vector3d normal{0., 0., 1.}, longitudinal{1., 0., 0.}, lateral{0., 1., 0.};
  double normal_scale{1.}, spin_lever_scale{1.};
  RoadContactFrame() = default;
  RoadContactFrame(double gradient_x, double gradient_y)
  {
    normal_scale = std::sqrt(1. + gradient_x * gradient_x + gradient_y * gradient_y);
    normal = {-gradient_x / normal_scale, -gradient_y / normal_scale, 1. / normal_scale};
    const double forward_norm = std::sqrt(1. + gradient_x * gradient_x);
    longitudinal = {1. / forward_norm, 0., gradient_x / forward_norm};
    lateral = normal.cross(longitudinal);
    spin_lever_scale = forward_norm / normal_scale;
  }
  Eigen::Vector2d velocity(const Eigen::Vector3d & hub_velocity) const
  {
    return {longitudinal.dot(hub_velocity), lateral.dot(hub_velocity)};
  }
  Eigen::Vector3d force(double fx, double fy) const { return fx * longitudinal + fy * lateral; }
};
}  // namespace tam::ocd::vehicle_dynamics
