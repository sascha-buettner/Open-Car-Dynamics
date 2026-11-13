# Optional coupled tire contact

Supply `tire.contact.enabled` and all physical contact parameters explicitly.
The public defaults are dummy placeholders and do not describe a vehicle. The enabled path uses five unilateral macro nodes
with a positive definite condensed elastic matrix (area diagonal plus graph
Laplacian). Precomputed active-set responses avoid a runtime FEM solve. This is
an engineering surrogate, not an identified rubber material or an ACC/iRacing
implementation. Width, footprint, coupling and contact damping require rig data.

Each node receives an enveloped height, gradient and heading derivative.
Normal force and yaw moment follow the contact potential. Tangential Magic
Formula forces use the local 3D tangent plane and a power-consistent wheel spin
lever. Contact acts on the unsprung system. The suspension transfers force to
the sprung body; road forces are not added directly to the chassis as well.

The optional vertical belt splits mass from the existing total unsprung mass;
it never adds that mass a second time. Belt/carcass and residual contact stiffness
preserve the requested total static tire stiffness. This vertical mode does not
model lateral or torsional ring dynamics. Full wheel spin inertia stays aggregated.
The four lateral relaxation states shift the patch using measured lateral stiffness;
the bounded sidewall height correction is a hypothesis, not a fitted tire law.

Optional nonlinear front/rear damper curves use wheel-to-damper motion ratio for
both velocity and reciprocal force. Curves must be finite, monotone, nonnegative
and start at the origin; extrapolation is logged. Optional progressive per-wheel
bump stops and third heave rubber transfer equal and opposite forces. Third
rubber preload is tared to preserve the nominal static equilibrium.

`strict_surface_coverage` rejects missing contact instead of silently reverting
to the ribbon. Configuration, missing/bridged nodes, normal/tangent forces,
contact energy, belt energy, deflection-range and damper extrapolation diagnostics
are logged. Tests cover complementarity, local energy gradients, contact-frame
power, belt mass/energy, heave energy and nonlinear damper passivity. Those local
checks do not certify a complete reciprocal lateral carcass material model.
