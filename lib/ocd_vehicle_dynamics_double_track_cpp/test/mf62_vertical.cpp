#include <cmath>
#include <iostream>
#include <stdexcept>

#include "ocd_vehicle_dynamics_double_track_cpp/mf62_vertical.hpp"
using namespace tam::ocd::vehicle_dynamics;

int main()
{
  // Synthetic tire: numbers define a valid numerical domain, not a vehicle fit.
  MF62Vertical p;
  p.UNLOADED_RADIUS = .30;
  p.RIM_RADIUS = .20;
  p.WIDTH = .20;
  p.ASPECT_RATIO = .50;
  p.LONGVL = 30.;
  p.FNOMIN = 3000.;
  p.NOMPRES = 100000.;
  p.VERTICAL_STIFFNESS = 200000.;
  p.Q_RE0 = 1.;
  p.Q_V1 = .001;
  p.Q_V2 = 0.;
  p.Q_FZ2 = .1;
  p.Q_FCX = p.Q_FCY = p.Q_FYS1 = p.Q_FYS2 = p.Q_FYS3 = 0.;
  p.Q_CAM1 = p.Q_CAM2 = p.Q_CAM3 = 0.;
  p.PFZ1 = .05;
  p.BREFF = .1;
  p.DREFF = .1;
  p.FREFF = .1;
  p.BOTTOM_OFFST = p.BOTTOM_STIFF = 0.;
  p.QSY1 = .01;
  p.QSY2 = p.QSY3 = p.QSY4 = p.QSY5 = p.QSY6 = p.QSY8 = 0.;
  p.QSY7 = 1.;
  for (double gauge : {90000., 100000., 110000.})
    for (double omega : {0., 100., 200.})
      for (double load : {1000., 3000., 5000.}) {
        const double radius = mf62_loaded_radius(p, load, gauge, omega, 0.);
        const double recovered = mf62_vertical_force(p, gauge, omega, 0., radius);
        if (!(std::isfinite(radius) && radius > p.RIM_RADIUS &&
              std::abs(recovered - load) < 1e-7 * load))
          throw std::runtime_error("MF6.2 radial inverse failed");
        if (!(mf62_effective_radius(p, load, gauge, omega) > p.RIM_RADIUS &&
              mf62_rolling_torque(p, load, gauge, omega * p.UNLOADED_RADIUS, 0., 0.) >= 0.))
          throw std::runtime_error("MF6.2 radius or rolling torque invalid");
      }
  if (!(mf62_free_radius(p, 200.) > mf62_free_radius(p, 100.) &&
        mf62_free_radius(p, 100.) > mf62_free_radius(p, 0.)))
    throw std::runtime_error("Free radius does not grow with speed");
  if (mf62_vertical_force(p, p.NOMPRES, 0., 0., p.UNLOADED_RADIUS + .01) != 0.)
    throw std::runtime_error("Airborne tire transmits normal force");
  std::cout << "Synthetic MF6.2 radial inverse, rolling torque and lift-off passed\n";
}
