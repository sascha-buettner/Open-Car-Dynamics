#pragma once

#include <array>
#include <string>
#include <vector>

// Public numerical fixtures only. These values describe no measured vehicle.
namespace synthetic_vehicle
{
template <class Model>
void contact(Model & model, bool belt = false)
{
  const auto parameters = model.get_param_manager();
  const std::string base = "vehicle_dynamics_double_track.";
  const auto set = [&](const std::string & name, double value) {
    parameters->set_value(base + name, value);
  };
  set("tire.contact.enabled", 1.);
  set("tire.contact.sampling_points", 5.);
  set("tire.contact.strict_surface_coverage", 0.);
  set("tire.contact.damping_Nspm.front", 0.);
  set("tire.contact.damping_Nspm.rear", 0.);
  set("tire.contact.width_m.front", .25);
  set("tire.contact.width_m.rear", .25);
  set("tire.contact.half_length_m", .08);
  set("tire.contact.enveloping_radius_m.front", .05);
  set("tire.contact.enveloping_radius_m.rear", .05);
  set("tire.contact.micro_radius_m", .01);
  set("tire.contact.lateral_stiffness_Npm.front", 100000.);
  set("tire.contact.lateral_stiffness_Npm.rear", 100000.);
  set("tire.contact.sidewall_height_m.front", .10);
  set("tire.contact.sidewall_height_m.rear", .10);
  set("tire.contact.lateral_drop_scale", 0.);
  set("tire.contact.belt.enabled", belt ? 1. : 0.);
  if (belt) {
    set("tire.contact.belt.mass_kg.front", 3.);
    set("tire.contact.belt.mass_kg.rear", 3.);
    set("tire.contact.belt.frequency_hz.front", 100.);
    set("tire.contact.belt.frequency_hz.rear", 100.);
    set("tire.contact.belt.damping_ratio.front", .1);
    set("tire.contact.belt.damping_ratio.rear", .1);
  }
}

template <class Model>
void geometry(Model & model)
{
  const auto parameters = model.get_param_manager();
  const std::string base = "vehicle_dynamics_double_track.suspension.geometry.";
  parameters->set_value(base + "travel_m", std::vector<double>{-.1, 0., .1});
  for (const char * axle : {"front", "rear"}) {
    parameters->set_value(base + "equivalent_travel_m." + axle,
                          std::vector<double>{-.1, 0., .1});
    parameters->set_value(base + "camber_delta_rad." + axle, std::vector<double>{0., 0., 0.});
    parameters->set_value(base + "toe_delta_rad." + axle, std::vector<double>{0., 0., 0.});
    parameters->set_value(base + "roll_center_height_m." + axle,
                          std::vector<double>{.2, .2, .2});
  }
  parameters->set_value(base + "steering_mean_angle_rad",
                        std::vector<double>{-.2, 0., .2});
  parameters->set_value(base + "jacking_body_m", std::vector<double>{0., 0., 0.});
  parameters->set_value(base + "enabled", 1.);
}

template <class Model>
void kinematics(Model & model)
{
  const auto parameters = model.get_param_manager();
  const std::string base = "vehicle_dynamics_double_track.suspension.kinematics.";
  parameters->set_value(base + "travel_m", std::vector<double>{-.1, .1});
  parameters->set_value(base + "steering_mean_angle_rad",
                        std::vector<double>{-.2, 0., .2});
  parameters->set_value(base + "rack_m", std::vector<double>{-.02, 0., .02});
  parameters->set_value("vehicle_dynamics_double_track.suspension.heave.rear_bar_travel_m",
                        std::vector<double>{-.1, .1});
  parameters->set_value("vehicle_dynamics_double_track.suspension.heave.rear_device_compression_m",
                        std::vector<double>{-.1, .1});
  const std::vector<double> identity{-.1, .1, -.1, .1};
  const auto constant = [](double v) { return std::vector<double>(4, v); };
  for (const char * axle : {"front", "rear"}) {
    const std::string suffix = std::string(".") + axle;
    parameters->set_value(base + "control_m" + suffix, std::vector<double>{-.02, .02});
    for (const char * name : {"spring", "damper", "arb", "rocker"})
      parameters->set_value(base + name + suffix, identity);
    parameters->set_value(base + "normal_x" + suffix, constant(0.));
    parameters->set_value(base + "normal_y" + suffix, constant(1.));
    parameters->set_value(base + "normal_z" + suffix, constant(0.));
    parameters->set_value(base + "contact_x" + suffix, constant(0.));
    parameters->set_value(base + "contact_y" + suffix, constant(0.));
    parameters->set_value(base + "contact_z" + suffix, constant(-.3));
    parameters->set_value(base + "jac_x" + suffix, constant(.1));
    parameters->set_value(base + "jac_y" + suffix, constant(.05));
    parameters->set_value(base + "jac_z" + suffix, constant(1.));
    for (const char * name : {"angjac_x", "angjac_y", "angjac_z"})
      parameters->set_value(base + name + suffix, constant(0.));
  }
  parameters->set_value(base + "enabled", 1.);
}
}  // namespace synthetic_vehicle
