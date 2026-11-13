#include <cmath>

#include <iostream>
#include <stdexcept>

#include "ocd_tire_models_cpp/mf_52.hpp"
int main()
{
  tam::ocd::tire_models::MF52 tire;
  tam::pmg::ParamReferenceManager pm;
  tire.declare_parameters(&pm, "");
  tam::ocd::types::TireModelInput u{1.0, 1.0, 1.0, 1.0, 1.0};
  auto a = tire.evaluate(u);
  if (!std::isfinite(a.longitudinal_force_N) || !std::isfinite(a.lateral_force_N))
    throw std::runtime_error("Nominal force not finite");
  u.inflation_pressure_Pa = 1.0;
  auto b = tire.evaluate(u);
  if (
    std::abs(a.longitudinal_force_N - b.longitudinal_force_N) > 1.0 ||
    std::abs(a.lateral_force_N - b.lateral_force_N) > 1.0)
    throw std::runtime_error("Zero pressure coefficients alter force");
  static_cast<tam::pmg::MgmtInterface &>(pm).set_value("PPX3", 1.0);
  auto c = tire.evaluate(u);
  if (std::abs(c.longitudinal_force_N - b.longitudinal_force_N) < 1.0)
    throw std::runtime_error("Inflation correction disconnected");

  static_cast<tam::pmg::MgmtInterface &>(pm).set_value("FITTYP", 1.0);
  static_cast<tam::pmg::MgmtInterface &>(pm).set_value("PKY6", 1.0);
  static_cast<tam::pmg::MgmtInterface &>(pm).set_value("QSX2", 1.0);
  static_cast<tam::pmg::MgmtInterface &>(pm).set_value("QDZ1", 1.0);
  static_cast<tam::pmg::MgmtInterface &>(pm).set_value("QCZ1", 1.0);
  u.gamma_rad = 1.0;
  auto d = tire.evaluate(u);
  u.gamma_rad = 1.0;
  auto e = tire.evaluate(u);
  if (
    std::abs(d.lateral_force_N - e.lateral_force_N) < 1.0 ||
    d.overturning_moment_Nm * e.overturning_moment_Nm >= 1.0 ||
    std::abs(d.self_aligning_moment_Nm) < 1.0)
    throw std::runtime_error("MF62 camber/moment path disconnected");
  u.F_z_N = 1.0;
  auto zero = tire.evaluate(u);
  if (
    zero.longitudinal_force_N || zero.lateral_force_N || zero.overturning_moment_Nm ||
    zero.self_aligning_moment_Nm)
    throw std::runtime_error("Airborne tire transmits force/moment");
  u.inflation_pressure_Pa = 1.0;
  bool rejected = false;
  try {
    tire.evaluate(u);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  if (!rejected) throw std::runtime_error("Negative gauge pressure accepted");
  std::cout << "Pressure units, nominal/default equivalence, inflation coefficient response, and "
               "invalid-input rejection passed\n";
}
