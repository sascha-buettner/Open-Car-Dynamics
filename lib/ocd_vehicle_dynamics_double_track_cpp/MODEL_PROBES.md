# Optional model probes

The synthetic C++ probes in the four OCD `lib/*/test/` directories are source-level development tools. They use generic, unmeasured inputs and are not a vehicle calibration or an Imola replay. Ordinary CMake and ROS builds do not compile them, even when `BUILD_TESTING=ON`.

From a sourced ROS workspace, build them only when needed:

```bash
colcon build --packages-select ocd_aerodynamics_models_cpp ocd_steering_actuator_pt1_cpp ocd_tire_models_cpp ocd_vehicle_dynamics_double_track_cpp --cmake-args -DOCD_BUILD_MODEL_PROBES=ON
```

Run the resulting executables directly or use CTest in those package build directories. The `OCD_BUILD_MODEL_PROBES` option defaults to `OFF`; turn it off again for ordinary builds.

| Probe group | Purpose |
|---|---|
| `positive_residual_root`, `hybrid_radial_root`, `coupled_radial_root`, `radial_root_contract`, `generic_radial_lookup` | Inspect radial root finding and contact stiffness near unload/recontact. |
| `longitudinal_relaxation`, `tire_load_memory`, `mf62_vertical`, `belt_energy` | Inspect the relaxation, load-memory and vertical-belt mechanisms implicated in braking recovery. |
| `five_point_contact`, `contact_energy`, `contact_energy_envelope`, `point_road_energy`, `contact_frame` | Check contact geometry and force/energy response on synthetic surfaces. |
| `suspension_geometry`, `coupled_geometry`, `heave_energy`, `heave_coil_energy`, `damper_curves` | Check suspension geometry and vertical force paths. |
| `plain_mf52_baseline`, `suspension_output`, `material_velocity`, `external_wheel_forces`, `mf_schwift_ring`, `mf_schwift_vehicle` | Check standalone and other model interfaces. |

The remaining tire, aero and steering probes live beside their owning packages and use the same opt-in flag. For real positive-slip episodes, use the private IAC Sim trajectory-replay diagnostic with an explicit IAC Launch profile and recorded inputs.
