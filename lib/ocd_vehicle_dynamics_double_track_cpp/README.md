# Double-track vehicle dynamics

The public model is unchanged in shape. Physical equations are separated under
`include/ocd_vehicle_dynamics_double_track_cpp/detail` into kinematics, radial
contact, tangential forces, suspension, body/aero, pipeline, derivatives and
read-only diagnostics. The detail headers are implementation units, not separate
stateful solvers. Each ODE evaluation uses the current stage state/forces.

Physical vehicle and tire profiles must be supplied by private configuration.
Public defaults and fixture values are placeholders. Independent
switches cover static front/rear roll centres, nonlinear suspension and coupled
kinematics, aero, radial law, force couplings, memory, belt and contact sampling.

One/five-point contact uses the same constitutive and structural tire. The
nonlinear radial-contact adapter extends the supplied MF6.2 radial law to either
sampling count without introducing a belt. Five-point samples remain massless.

The rolling radial lookup preserves the distinction between tangent stiffness
and its integrated compliance/secant force coefficient. The safeguarded TOMS748
fallback solves the original positive constitutive relation; belt damping uses
the consistent tangent at each trial. The original force-domain guards remain.
Smaller integration steps cannot repair an invalid constitutive law.

`get_stationary_diagnostics()` is a read-only typed measurement surface. The
historical packed array encoding belongs only to the evaluation bench. Native
force/energy, geometry, root-solver and ring-kernel contracts live in `test/`.

`mf_schwift_ring.hpp` is an independently implemented six-coordinate linearized
ring kernel, not an alias of the vertical belt and not a complete vehicle-admitted
vendor MF-Swift implementation. Full vehicle/contact/spin/inertia coupling has
an explicit separate admission gate. Do not call the steady MF6.2 branch Swift.
