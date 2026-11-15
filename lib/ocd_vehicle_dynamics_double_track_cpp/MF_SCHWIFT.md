# mfSCHWIFT: independent open rigid-ring tire

The `tire.mfschwift.enabled` switch selects a vehicle-coupled, small-displacement
six-coordinate rigid ring. It is an independent implementation, not vendor
MF-Swift code or an assertion of proprietary model equivalence. The existing
vertical-belt model remains a separate branch; steady MF6.2 alone is not a ring.

## Coordinates and force paths

- Each wheel has ring translations x/y/z, camber, windup and yaw. Vertical
  translation reuses the mass-partitioned belt state; the other five coordinates
  and their rates add ten states per wheel.
- x/y displacement is relative to the hub in the moving body frame. Planar ring
  velocity is absolute, expressed in that frame. Modal angular velocities are
  absolute in nonrotating hub axes; windup rate subtracts hub angular velocity
  and rotor speed. This avoids differentiating delayed spindle commands.
- Contact friction and the horizontal reactions of scanned-road normals act on
  the ring. Only sidewall reactions reach the spindle/body; no second direct
  contact-force path is added. Vertical contact reaches the ring/belt and then
  the wheel mass through the radial sidewall.
- A circular ring is invariant under spin. Conservative normal enveloping has
  no windup-coordinate torque. Tangential traction and MF rolling resistance
  act on spin. Camber rotates the patch about the ring centre, moving its query
  centre laterally by R times camber and its node heights by camber times y.
  Sample velocities include relative ring translation, camber and yaw motion.
- The body planar mass matrix excludes ring translational mass and reference
  yaw inertia. Ring transverse inertia is removed from body roll/pitch inertia.
  The iac_sim profile removes ring spin inertia from the wheel rotor inertia.
  Partition guards reject nonpositive masses or inertia matrices.

## Structural and residual compliance

Diagonal support stiffnesses use k = m(2*pi*f)^2 or I(2*pi*f)^2; damping uses
c = 2*zeta*m*2*pi*f or its rotational counterpart. The coupled lateral energy is

    U_lat = Ky*(y + a*R*camber)^2/2 + Kcamber*camber^2/2.

The support forces and moments are derivatives of this same energy, with a
matching positive dissipative quadratic form. The shear-centre location and
modal frequencies must be supplied by private configuration.

The quasi-static carcass compliance at the ground is

    Jx = 1/Kx + R^2/Kwindup
    Jy = 1/Ky + ((1-a)*R)^2/Kcamber.

For supplied total tangential stiffness Ktotal, residual compliance is
1/Kres = 1/Ktotal - J. Nonpositive residual compliance is rejected. Local native
MF6.2 slopes Ckappa/Calpha give residual relaxation time C/(Kres*|Vx|).
The speed, load, and time floors are configurable extrapolation assumptions.
They are not identified racing-tire dynamics at near-zero load or speed.

The gyroscopic pair is skew and uses the ring spin inertia and signed spin.
Rolling resistance reverses smoothly at zero spin instead of applying a finite
constant torque to a stationary ring.

## Physical parameters

All vehicle and tire values, including ring mass, inertia, modal frequencies,
damping, stiffness, pressure, and radial lookup data, must be supplied through
private configuration. Public defaults are deliberately dummy placeholders.

## Readouts and limits

The model is an open small-displacement approximation. It does not establish exact
nonlinear global passivity or measured tire calibration.

`kinematic/wheel_longitudinal_slip` is the hub-wheel encoder-equivalent physical
slip. `imr/tire_longitudinal_slip` is ring/contact kinematic slip when this branch
is active; `kappa_tire_*` is its relaxed force state. Do not compare the latter two
to measured wheel slip without identifying the different quantity. Controller
CAN wheel-speed routing is unchanged.

## Equation coverage and limits

This is a small-angle modal ring with nonrotating support coordinates, positive
residual compliance, MF6.2 steady forces/moments and spatial normal contact.
It does not implement proprietary rotating-damper corrections, finite-angle
SO(3) transport, flexible-shell modes, exact proprietary ellipsoidal cams or a
proprietary residual dynamic kernel. Supplied structural values do not validate
these approximations or uniquely identify damping and contact dynamics.

Primary structural references:

- [op het Veld2008, Enhancing the MF-Swift Tyre Model for Inflation Pressure Changes](https://pure.tue.nl/ws/files/46975582/674114-1.pdf), printed8/PDF21, equations2.13–2.18.
- [Besselink, Schmeitz and Pacejka2010](https://doi.org/10.1080/00423111003748088), structural/contact compliance decomposition.
- [MBDyn1.7.3 input manual](https://www.mbdyn.org/userfiles/documents/mbdyn-input-1.7.3.pdf), wheel4 open ring/patch interface reference. No implementation from MBDyn is copied here.


## Safeguarded radial algebraic closure

The original converged force/coefficient iteration remains the first path. Its
fallback solves the original positive force-coefficient equation with TOMS748.
If tangent feedback is not contractive, a second positive scalar solve closes the
original tangent equation. The effective-radius/contact-geometry feedback also
retains its converged fixed point, with a radius-root fallback bounded above by
the unloaded MF radius at the current spin. Each nested trial resets its radius
seed; the final original residual and unchanged force/domain guards are checked.
A narrow bracket across a discontinuity is rejected when its residual is invalid.
These changes do not clip the force law or prove a root exists in every state.

MF steady camber coefficients currently receive the configured spindle/kinematic
camber. Ring camber enters the explicit contact geometry and modal motion; a
separately identified residual-camber force model is outside the validated scope.


When the nested scalar fallbacks cannot close all three equations, the v8
fallback solves secant stiffness, tangent stiffness and effective radius
simultaneously with bounded Newton/Levenberg steps. The v9 fallback then tries
Eigen's MINPACK-style Powell hybrid in positive bounded logistic coordinates.
Two search scalings are allowed, but acceptance always checks the original
coupled residual (in its declared normalized coordinates) below 3e-10. A small
step or a solver termination status alone is not a certificate. Existing
positive-force and road-coverage guards remain active. A failed search does not
prove that a physical root is absent.

Driving results and their exact parameters are retained in the private
evaluation repository. This public document does not contain run data.

- [Alefeld, Potra and Shi 1995](https://doi.org/10.1145/210089.210111), safeguarded scalar TOMS748 bracketing.
- [Eigen nonlinear optimization module](https://libeigen.gitlab.io/eigen/docs-nightly/unsupported/group__NonLinearOptimization__Module.html), MINPACK port and Powell hybrid.
- [Original MINPACK HYBRD manual](https://www.math.utah.edu/software/minpack/minpack/hybrd.html), trust-region root search and termination limits.
