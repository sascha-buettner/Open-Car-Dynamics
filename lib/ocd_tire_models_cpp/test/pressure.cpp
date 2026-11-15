#include <cmath>
#include <stdexcept>

#include "ocd_tire_models_cpp/mf_52.hpp"

int main()
{
  tam::ocd::tire_models::MF52 tire;
  tam::pmg::ParamReferenceManager parameters;
  tire.declare_parameters(&parameters, "");
  tam::ocd::types::TireModelInput input{.10, .04, 3000., .02, 1.};
  const auto nominal = tire.evaluate(input);
  if (!std::isfinite(nominal.longitudinal_force_N) ||
      !std::isfinite(nominal.lateral_force_N))
    throw std::runtime_error("Public MF5.2 baseline is not finite");
  input.inflation_pressure_Pa = 1.2;
  const auto neutral = tire.evaluate(input);
  if (std::abs(neutral.longitudinal_force_N - nominal.longitudinal_force_N) > 1e-8 ||
      std::abs(neutral.lateral_force_N - nominal.lateral_force_N) > 1e-8)
    throw std::runtime_error("Unconfigured inflation correction changed MF5.2");
  static_cast<tam::pmg::MgmtInterface &>(parameters).set_value("PPX3", .5);
  const auto pressure_sensitive = tire.evaluate(input);
  if (std::abs(pressure_sensitive.longitudinal_force_N - nominal.longitudinal_force_N) < 1.)
    throw std::runtime_error("Opt-in inflation correction is disconnected");

  static_cast<tam::pmg::MgmtInterface &>(parameters).set_value("FITTYP", 62.);
  static_cast<tam::pmg::MgmtInterface &>(parameters).set_value("UNLOADED_RADIUS", .3);
  static_cast<tam::pmg::MgmtInterface &>(parameters).set_value("FZMIN", 100.);
  const auto mf62 = tire.evaluate(input);
  if (!std::isfinite(mf62.longitudinal_force_N) ||
      !std::isfinite(mf62.lateral_force_N) ||
      !std::isfinite(mf62.self_aligning_moment_Nm) ||
      !std::isfinite(mf62.overturning_moment_Nm))
    throw std::runtime_error("Opt-in MF6.2 path is not finite");
  input.F_z_N = 0.;
  const auto airborne = tire.evaluate(input);
  if (airborne.longitudinal_force_N != 0. || airborne.lateral_force_N != 0. ||
      airborne.overturning_moment_Nm != 0. || airborne.self_aligning_moment_Nm != 0.)
    throw std::runtime_error("Airborne MF6.2 tire transmits force");
  input.inflation_pressure_Pa = 0.;
  bool rejected = false;
  try {
    tire.evaluate(input);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Nonpositive gauge pressure accepted");
}
