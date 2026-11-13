// Copyright 2026 Simon Sagmeister
#pragma once

#include <algorithm>
#include <boost/math/tools/toms748_solve.hpp>
#include <limits>
#include <string>

#include "ocd_vehicle_dynamics_double_track_cpp/stage_road_cache.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/vehicle_dynamics_equations.hpp"

#define EVAL_MACRO_PER_WHEEL(x) \
  x(front_left);                \
  x(front_right);               \
  x(rear_left);                 \
  x(rear_right)

// The implementation is grouped by physical subsystem. Public parameters and
// state ordering remain defined by the model interface.
#include "ocd_vehicle_dynamics_double_track_cpp/detail/body_dynamics.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/detail/diagnostics.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/detail/equation_pipeline.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/detail/kinematics.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/detail/radial_contact.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/detail/state_derivatives.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/detail/suspension.hpp"
#include "ocd_vehicle_dynamics_double_track_cpp/detail/tangential_forces.hpp"

#undef EVAL_MACRO_PER_WHEEL
