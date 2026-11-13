// Copyright 2026 Simon Sagmeister

#include "ocd_tire_models_cpp/mf_52.hpp"

#include <algorithm>
#include <stdexcept>
namespace tam::ocd::tire_models
{
types::TireModelOutput MF52::evaluate(types::TireModelInput const & input)
{
  if (p_.FITTYP == 62.) return evaluate62(input);
  types::TireModelOutput out;
  if (!(std::isfinite(input.inflation_pressure_Pa) && input.inflation_pressure_Pa > 0. &&
        p_.NOMPRES > 0.))
    throw std::invalid_argument("Tire gauge pressure and nominal pressure must be positive");
  const double dp = (input.inflation_pressure_Pa - p_.NOMPRES) / p_.NOMPRES;
  double F_x_o = long_force_pure_slip(input.F_z_N, input.longitudinal_slip, dp);
  double F_y_o = lat_force_pure_slip(input.F_z_N, -input.slip_angle_rad, input.gamma_rad, dp);
  out.longitudinal_force_N = long_force_combined_slip(
    input.F_z_N, input.longitudinal_slip, -input.slip_angle_rad, F_x_o, input.gamma_rad);
  out.lateral_force_N = lat_force_combined_slip(
    input.F_z_N, input.longitudinal_slip, -input.slip_angle_rad, F_y_o, input.gamma_rad);
  // Implement MF52 in here
  return out;
}
double MF52::long_force_pure_slip(double F_z, double kappa, double dp) const
{
  double F_z_o_a = p_.LFZO * p_.FNOMIN;
  double dfz = (F_z - F_z_o_a) / F_z_o_a;
  double S_H_x = (p_.PHX1 + p_.PHX2 * dfz) * p_.LHX;
  double S_V_x = F_z * (p_.PVX1 + p_.PVX2 * dfz) * p_.LVX * p_.LMUX;
  double kappa_x = kappa + S_H_x;
  double C_x = p_.PCX1 * p_.LCX;
  double mu_x = (p_.PDX1 + p_.PDX2 * dfz) * (1. + p_.PPX3 * dp + p_.PPX4 * dp * dp) * p_.LMUX;
  double D_x = mu_x * F_z;
  double E_x = (p_.PEX1 + p_.PEX2 * dfz + p_.PEX3 * std::pow(dfz, 2)) *
               (1.0 - p_.PEX4 * std::copysign(1.0, kappa_x)) * p_.LEX;
  double K_x = F_z * (p_.PKX1 + p_.PKX2 * dfz) * std::exp(p_.PKX3 * dfz) *
               (1. + p_.PPX1 * dp + p_.PPX2 * dp * dp) * p_.LKX;
  double B_x = K_x / (C_x * D_x);
  // Longitudinal Force (pure longitudinal slip)
  double F_x_o =
    D_x *
      std::sin(C_x * std::atan(B_x * kappa_x - E_x * (B_x * kappa_x - std::atan(B_x * kappa_x)))) +
    S_V_x;
  return F_x_o;
}
double MF52::lat_force_pure_slip(double F_z, double alpha, double gamma, double dp) const
{
  double gamma_s = std::sin(gamma);
  double F_z_o_a = p_.LFZO * p_.FNOMIN;
  double dfz = (F_z - F_z_o_a) / F_z_o_a;
  double K_y_g_o = F_z * (p_.PKY6 + p_.PKY7 * dfz) * (1. + p_.PPY5 * dp) * p_.LKYG;
  double S_V_y_g = F_z * (p_.PVY3 + p_.PVY4 * dfz) * gamma_s * p_.LKYG * p_.LMUY;
  double S_V_y = F_z * (p_.PVY1 + p_.PVY2 * dfz) * p_.LVY * p_.LMUY + S_V_y_g;
  double K_y =
    p_.PKY1 * F_z_o_a * (1.0 - p_.PKY3 * std::abs(gamma_s)) *
    std::sin(
      p_.PKY4 *
      std::atan(
        F_z / ((p_.PKY2 + p_.PKY5 * std::pow(gamma_s, 2)) * F_z_o_a * (1. + p_.PPY2 * dp)))) *
    (1. + p_.PPY1 * dp) * p_.LKY;
  double S_H_y = (p_.PHY1 + p_.PHY2 * dfz) * p_.LHY + (K_y_g_o * gamma_s - S_V_y_g) / K_y;
  double C_y = p_.PCY1 * p_.LCY;
  double mu_y = ((p_.PDY1 + p_.PDY2 * dfz) * (1.0 - p_.PDY3 * std::pow(gamma_s, 2))) *
                (1. + p_.PPY3 * dp + p_.PPY4 * dp * dp) * p_.LMUY;
  double D_y = mu_y * F_z;
  double B_y = K_y / (C_y * D_y);
  double alpha_y = alpha + S_H_y;
  double E_y = (p_.PEY1 + p_.PEY2 * dfz) *
               (1.0 + p_.PEY5 * std::pow(gamma_s, 2) -
                (p_.PEY3 + p_.PEY4 * gamma_s) * std::copysign(1.0, alpha_y)) *
               p_.LEY;
  double F_y_o =
    D_y *
      std::sin(C_y * std::atan(B_y * alpha_y - E_y * (B_y * alpha_y - std::atan(B_y * alpha_y)))) +
    S_V_y;
  return F_y_o;
}
double MF52::long_force_combined_slip(
  double F_z, double kappa, double alpha, double F_x_o, double gamma) const
{
  double gamma_s = std::sin(gamma);
  double F_z_o_a = p_.LFZO * p_.FNOMIN;
  double dfz = (F_z - F_z_o_a) / F_z_o_a;
  double S_H_x_a = p_.RHX1;
  double E_x_a = p_.REX1 + p_.REX2 * dfz;
  double C_x_a = p_.RCX1;
  double B_x_a =
    (p_.RBX1 + p_.RBX3 * std::pow(gamma_s, 2)) * std::cos(std::atan(p_.RBX2 * kappa)) * p_.LXAL;
  double alpha_S = alpha + S_H_x_a;
  double G_x_a_o = std::cos(
    C_x_a * std::atan(B_x_a * S_H_x_a - E_x_a * (B_x_a * S_H_x_a - std::atan(B_x_a * S_H_x_a))));
  double G_x_a =
    std::cos(
      C_x_a * std::atan(B_x_a * alpha_S - E_x_a * (B_x_a * alpha_S - std::atan(B_x_a * alpha_S)))) /
    G_x_a_o;
  double F_x = G_x_a * F_x_o;
  return F_x;
}
double MF52::lat_force_combined_slip(
  double F_z, double kappa, double alpha, double F_y_o, double gamma) const
{
  double gamma_s = std::sin(gamma);
  double F_z_o_a = p_.LFZO * p_.FNOMIN;
  double dfz = (F_z - F_z_o_a) / F_z_o_a;
  double mu_y = ((p_.PDY1 + p_.PDY2 * dfz) / (1.0 + p_.PDY3 * std::pow(gamma_s, 2))) * p_.LMUY;
  double D_V_y_k = mu_y * F_z * (p_.RVY1 + p_.RVY2 * dfz + p_.RVY3 * gamma_s) *
                   std::cos(std::atan(p_.RVY4 * alpha));
  double S_V_y_k = D_V_y_k * std::sin(p_.RVY5 * std::atan(p_.RVY6 * kappa)) * p_.LVYKA;
  double S_H_y_k = p_.RHY1 + p_.RHY2 * dfz;
  double E_y_k = p_.REY1 + p_.REY2 * dfz;
  double C_y_k = p_.RCY1;
  double B_y_k = (p_.RBY1 + p_.RBY4 * std::pow(gamma_s, 2)) *
                 std::cos(std::atan(p_.RBY2 * (alpha - p_.RBY3))) * p_.LYKA;
  double kappa_S = kappa + S_H_y_k;
  double G_y_k_o = std::cos(
    C_y_k * std::atan(B_y_k * S_H_y_k - E_y_k * (B_y_k * S_H_y_k - std::atan(B_y_k * S_H_y_k))));
  double G_y_k =
    std::cos(
      C_y_k * std::atan(B_y_k * kappa_S - E_y_k * (B_y_k * kappa_S - std::atan(B_y_k * kappa_S)))) /
    G_y_k_o;
  double F_y = G_y_k * F_y_o + S_V_y_k;
  return F_y;
}
void MF52::declare_parameters(
  tam::pmg::ParamReferenceManager * param_manager, std::string name_prefix)
{
  param_manager->declare_parameter(
    name_prefix + "FITTYP", &p_.FITTYP, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "UNLOADED_RADIUS", &p_.UNLOADED_RADIUS, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "FZMIN", &p_.FZMIN, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "LKYC", &p_.LKYC, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "LKZC", &p_.LKZC, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX1", &p_.QSX1, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX2", &p_.QSX2, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX3", &p_.QSX3, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX4", &p_.QSX4, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX5", &p_.QSX5, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX6", &p_.QSX6, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX7", &p_.QSX7, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX8", &p_.QSX8, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX9", &p_.QSX9, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX10", &p_.QSX10, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX11", &p_.QSX11, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX12", &p_.QSX12, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX13", &p_.QSX13, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QSX14", &p_.QSX14, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "PPMX1", &p_.PPMX1, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QBZ1", &p_.QBZ1, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QBZ2", &p_.QBZ2, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QBZ3", &p_.QBZ3, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QBZ4", &p_.QBZ4, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QBZ5", &p_.QBZ5, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QBZ9", &p_.QBZ9, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QBZ10", &p_.QBZ10, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QCZ1", &p_.QCZ1, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QDZ1", &p_.QDZ1, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QDZ2", &p_.QDZ2, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QDZ3", &p_.QDZ3, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QDZ4", &p_.QDZ4, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QDZ6", &p_.QDZ6, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QDZ7", &p_.QDZ7, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QDZ8", &p_.QDZ8, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QDZ9", &p_.QDZ9, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QDZ10", &p_.QDZ10, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QDZ11", &p_.QDZ11, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QEZ1", &p_.QEZ1, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QEZ2", &p_.QEZ2, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QEZ3", &p_.QEZ3, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QEZ4", &p_.QEZ4, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QEZ5", &p_.QEZ5, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QHZ1", &p_.QHZ1, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QHZ2", &p_.QHZ2, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QHZ3", &p_.QHZ3, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "QHZ4", &p_.QHZ4, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "PPZ1", &p_.PPZ1, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "PPZ2", &p_.PPZ2, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "SSZ1", &p_.SSZ1, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "SSZ2", &p_.SSZ2, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "SSZ3", &p_.SSZ3, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  param_manager->declare_parameter(
    name_prefix + "SSZ4", &p_.SSZ4, 1.0, tam::pmg::ParameterType::DOUBLE,
    "MF6.2 manufacturer coefficient; opt-in FITTYP62");
  // clang-format off
  // Long parameters
  param_manager->declare_parameter(name_prefix + "PCX1", &(p_.PCX1), 1.6    , tam::pmg::ParameterType::DOUBLE, "Shape factor Cfx for longitudinal force"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PDX1", &(p_.PDX1), 1.5    , tam::pmg::ParameterType::DOUBLE, "Longitudinal friction Mux at Fznom"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PDX2", &(p_.PDX2), -0.04  , tam::pmg::ParameterType::DOUBLE, "Variation of friction Mux with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PDX3", &(p_.PDX3), 0.1    , tam::pmg::ParameterType::DOUBLE, "Variation of friction Mux with camber"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PEX1", &(p_.PEX1), 0.7    , tam::pmg::ParameterType::DOUBLE, "Longitudinal curvature Efx at Fznom"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PEX2", &(p_.PEX2), -0.17  , tam::pmg::ParameterType::DOUBLE, "Variation of curvature Efx with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PEX3", &(p_.PEX3), 0.023  , tam::pmg::ParameterType::DOUBLE, "Variation of curvature Efx with load squared"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PEX4", &(p_.PEX4), -0.14  , tam::pmg::ParameterType::DOUBLE, "Factor in curvature Efx while driving"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PKX1", &(p_.PKX1), 30.7   , tam::pmg::ParameterType::DOUBLE, "Longitudinal slip stiffness Kfx/Fz at Fznom"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PKX2", &(p_.PKX2), 0.27   , tam::pmg::ParameterType::DOUBLE, "Variation of slip stiffness Kfx/Fz with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PKX3", &(p_.PKX3), 0.13   , tam::pmg::ParameterType::DOUBLE, "Exponent in slip stiffness Kfx/Fz with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PHX1", &(p_.PHX1), 0.0    , tam::pmg::ParameterType::DOUBLE, "Horizontal shift Shx at Fznom"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PHX2", &(p_.PHX2), 0.0    , tam::pmg::ParameterType::DOUBLE, "Variation of shift Shx with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PVX1", &(p_.PVX1), 0.0    , tam::pmg::ParameterType::DOUBLE, "Vertical shift Svx/Fz at Fznom"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PVX2", &(p_.PVX2), 0.0    , tam::pmg::ParameterType::DOUBLE, "ariation of shift Svx/Fz with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RBX1", &(p_.RBX1), 17.4   , tam::pmg::ParameterType::DOUBLE, "Slope factor for combined slip Fx reduction"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RBX2", &(p_.RBX2), 12.9   , tam::pmg::ParameterType::DOUBLE, "Variation of slope Fx reduction with kappa"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RBX3", &(p_.RBX3), 0.0    , tam::pmg::ParameterType::DOUBLE, "Camber influence on combined slip longitudinal force"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RCX1", &(p_.RCX1), 1.1    , tam::pmg::ParameterType::DOUBLE, "Shape factor for combined slip Fx reduction"); // NOLINT
  param_manager->declare_parameter(name_prefix + "REX1", &(p_.REX1), 0.0    , tam::pmg::ParameterType::DOUBLE, "Curvature factor of combined Fx"); // NOLINT
  param_manager->declare_parameter(name_prefix + "REX2", &(p_.REX2), 0.0    , tam::pmg::ParameterType::DOUBLE, "Curvature factor of combined Fx with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RHX1", &(p_.RHX1), 0.001  , tam::pmg::ParameterType::DOUBLE, "Shift factor for combined slip Fx reduction"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PTX1", &(p_.PTX1), 0.75   , tam::pmg::ParameterType::DOUBLE, "Relaxation length SigKap0/Fz at Fznom"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PTX2", &(p_.PTX2), 0.007  , tam::pmg::ParameterType::DOUBLE, "Variation of SigKap0/Fz with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PTX3", &(p_.PTX3), -0.013 , tam::pmg::ParameterType::DOUBLE, "Variation of SigKap0/Fz with exponent of load"); // NOLINT
  // clang-format on

  // clang-format off
  // Lat parameters
  param_manager->declare_parameter(name_prefix + "PCY1", &(p_.PCY1), 1.5    , tam::pmg::ParameterType::DOUBLE, "Shape factor Cfy for lateral forces"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PDY1", &(p_.PDY1), 1.2    , tam::pmg::ParameterType::DOUBLE, "Lateral friction Muy "); // NOLINT
  param_manager->declare_parameter(name_prefix + "PDY2", &(p_.PDY2), -0.09  , tam::pmg::ParameterType::DOUBLE, "Variation of friction Muy with load "); // NOLINT
  param_manager->declare_parameter(name_prefix + "PDY3", &(p_.PDY3), 0.1    , tam::pmg::ParameterType::DOUBLE, "Variation of friction Muy with squared camber"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PEY1", &(p_.PEY1), 0.1    , tam::pmg::ParameterType::DOUBLE, "Lateral curvature Efy at Fznom"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PEY2", &(p_.PEY2), -1.2   , tam::pmg::ParameterType::DOUBLE, "Variation of curvature Efy with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PEY3", &(p_.PEY3), 0.05   , tam::pmg::ParameterType::DOUBLE, "Zero order camber dependency of curvature Efy"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PEY4", &(p_.PEY4), -8.1   , tam::pmg::ParameterType::DOUBLE, "Variation of curvature Efy with camber"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PKY1", &(p_.PKY1), -75.5  , tam::pmg::ParameterType::DOUBLE, "Maximum value of stiffness Kfy/Fznom"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PKY2", &(p_.PKY2), 4.65   , tam::pmg::ParameterType::DOUBLE, "Load at which Kfy reaches maximum value"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PKY3", &(p_.PKY3), -0.6   , tam::pmg::ParameterType::DOUBLE, "Variation of Kfy/Fznom with camber"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PHY1", &(p_.PHY1), 0.003  , tam::pmg::ParameterType::DOUBLE, "Horizontal shift Shy at Fznom"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PHY2", &(p_.PHY2), -0.0006, tam::pmg::ParameterType::DOUBLE, "Variation of shift Shy with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PHY3", &(p_.PHY3), 0.004  , tam::pmg::ParameterType::DOUBLE, "Variation of shift Shy with camber"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PVY1", &(p_.PVY1), 0.04   , tam::pmg::ParameterType::DOUBLE, "Vertical shift in Svy/Fz at Fznom"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PVY2", &(p_.PVY2), -0.02  , tam::pmg::ParameterType::DOUBLE, "Variation of shift Svy/Fz with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PVY3", &(p_.PVY3), -0.97  , tam::pmg::ParameterType::DOUBLE, "Variation of shift Svy/Fz with camber"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PVY4", &(p_.PVY4), -0.41  , tam::pmg::ParameterType::DOUBLE, "Variation of shift Svy/Fz with camber and load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RBY1", &(p_.RBY1), 20.6   , tam::pmg::ParameterType::DOUBLE, "Slope factor for combined Fy reduction"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RBY2", &(p_.RBY2), -23.5  , tam::pmg::ParameterType::DOUBLE, "Variation of slope Fy reduction with alpha"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RBY3", &(p_.RBY3), 0.001  , tam::pmg::ParameterType::DOUBLE, "Shift term for alpha in slope Fy reduction"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RCY1", &(p_.RCY1), 1.0    , tam::pmg::ParameterType::DOUBLE, "Shape factor for combined Fy reduction"); // NOLINT
  param_manager->declare_parameter(name_prefix + "REY1", &(p_.REY1), 0.0    , tam::pmg::ParameterType::DOUBLE, "Curvature factor of combined Fy"); // NOLINT
  param_manager->declare_parameter(name_prefix + "REY2", &(p_.REY2), 0.0    , tam::pmg::ParameterType::DOUBLE, "Curvature factor of combined Fy with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RHY1", &(p_.RHY1), -0.02  , tam::pmg::ParameterType::DOUBLE, "Shift factor for combined Fy reduction"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RHY2", &(p_.RHY2), 0.0    , tam::pmg::ParameterType::DOUBLE, "Shift factor for combined Fy reduction with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RVY1", &(p_.RVY1), 0.16   , tam::pmg::ParameterType::DOUBLE, "Kappa induced side force Svyk/Muy*Fz at Fznom"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RVY2", &(p_.RVY2), 0.03   , tam::pmg::ParameterType::DOUBLE, "Variation of Svyk/Muy*Fz with load"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RVY3", &(p_.RVY3), 27.5   , tam::pmg::ParameterType::DOUBLE, "Variation of Svyk/Muy*Fz with camber"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RVY4", &(p_.RVY4), -29.7  , tam::pmg::ParameterType::DOUBLE, "Variation of Svyk/Muy*Fz with alpha"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RVY5", &(p_.RVY5), 0.03   , tam::pmg::ParameterType::DOUBLE, "Variation of Svyk/Muy*Fz with kappa"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RVY6", &(p_.RVY6), 0.0    , tam::pmg::ParameterType::DOUBLE, "Variation of Svyk/Muy*Fz with atan(kappa)"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PTY1", &(p_.PTY1), 3.8    , tam::pmg::ParameterType::DOUBLE, "Peak value of relaxation length SigAlp0/R0"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PTY2", &(p_.PTY2), 4.7    , tam::pmg::ParameterType::DOUBLE, "Value of Fz/Fznom where SigAlp0 is extreme"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PEY5", &(p_.PEY5), 0.0    , tam::pmg::ParameterType::DOUBLE, "Lateral curvature variation with camber squared"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PKY4", &(p_.PKY4), 2.0    , tam::pmg::ParameterType::DOUBLE, "Lateral force stiffness, KFy curvature"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PKY5", &(p_.PKY5), 0.0    , tam::pmg::ParameterType::DOUBLE, "Variation of peak stiffness with squared camber"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PKY6", &(p_.PKY6), 0.0    , tam::pmg::ParameterType::DOUBLE, "Lateral force, Fy, camber stiffness factor"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PKY7", &(p_.PKY7), 0.0    , tam::pmg::ParameterType::DOUBLE, "Camber stiffness vertical load dependency"); // NOLINT
  param_manager->declare_parameter(name_prefix + "NOMPRES", &(p_.NOMPRES), 1.0, tam::pmg::ParameterType::DOUBLE, "MF6 inflation correction from supplied tire file");
  param_manager->declare_parameter(name_prefix + "PPX1", &(p_.PPX1), 1.0, tam::pmg::ParameterType::DOUBLE, "MF6 inflation correction from supplied tire file");
  param_manager->declare_parameter(name_prefix + "PPX2", &(p_.PPX2), 1.0, tam::pmg::ParameterType::DOUBLE, "MF6 inflation correction from supplied tire file");
  param_manager->declare_parameter(name_prefix + "PPX3", &(p_.PPX3), 1.0, tam::pmg::ParameterType::DOUBLE, "MF6 inflation correction from supplied tire file");
  param_manager->declare_parameter(name_prefix + "PPX4", &(p_.PPX4), 1.0, tam::pmg::ParameterType::DOUBLE, "MF6 inflation correction from supplied tire file");
  param_manager->declare_parameter(name_prefix + "PPY1", &(p_.PPY1), 0.0    , tam::pmg::ParameterType::DOUBLE, "Cornering stiffness variation with inflation pressure"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PPY2", &(p_.PPY2), 0.0    , tam::pmg::ParameterType::DOUBLE, "Cornering stiffness variation with inflation pressure induced nominal load dependency"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PPY3", &(p_.PPY3), 0.0    , tam::pmg::ParameterType::DOUBLE, "Linear inflation pressure on peak lateral friction"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PPY4", &(p_.PPY4), 0.0    , tam::pmg::ParameterType::DOUBLE, "Quadratic inflation pressure on peak lateral friction"); // NOLINT
  param_manager->declare_parameter(name_prefix + "PPY5", &(p_.PPY5), 0.0    , tam::pmg::ParameterType::DOUBLE, "Inflation pressure effect on camber stiffness"); // NOLINT
  param_manager->declare_parameter(name_prefix + "RBY4", &(p_.RBY4), 0.0    , tam::pmg::ParameterType::DOUBLE, "Lateral force, Fy, combined stiffness variation from camber"); // NOLINT
  // clang-format on

  // clang-format off
  // Scaling factors
  param_manager->declare_parameter(name_prefix +  "LFZO" , &(p_.LFZO)  , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of nominal (rated) load"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LCX"  , &(p_.LCX)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fx shape factor"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LMUX" , &(p_.LMUX)  , 0.97 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fx peak friction coefficient"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LEX"  , &(p_.LEX)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fx curvature factor"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LKX"  , &(p_.LKX)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fx slip stiffness"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LHX"  , &(p_.LHX)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fx horizontal shift"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LVX"  , &(p_.LVX)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fx vertical shift"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LGAX" , &(p_.LGAX)  , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of camber for Fx"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LCY"  , &(p_.LCY)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fy shape factor"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LMUY" , &(p_.LMUY)  , 0.97 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fy peak friction coefficient"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LEY"  , &(p_.LEY)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fy curvature factor"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LKY"  , &(p_.LKY)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fy cornering stiffness"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LKYG" , &(p_.LKYG)  , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fy camber force stiffness"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LHY"  , &(p_.LHY)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fy horizontal shift"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LVY"  , &(p_.LVY)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Fy vertical shift"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LGAY" , &(p_.LGAY)  , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of camber for Fy"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LTR"  , &(p_.LTR)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Peak of pneumatic trail"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LRES" , &(p_.LRES)  , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor for offset of residual torque"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LGAZ" , &(p_.LGAZ)  , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of camber for Mz"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LXAL" , &(p_.LXAL)  , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of alpha influence on Fx"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LYKA" , &(p_.LYKA)  , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of alpha influence on Fx"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LVYKA", &(p_.LVYKA) , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of kappa induced Fy"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LS"   , &(p_.LS)    , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Moment arm of Fx"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LSGKP", &(p_.LSGKP) , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Relaxation length of Fx"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LSGAL", &(p_.LSGAL) , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Relaxation length of Fy"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LGYR" , &(p_.LGYR)  , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of gyroscopic torque"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LMX"  , &(p_.LMX)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of overturning couple"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LVMX" , &(p_.LVMX)  , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of Mx vertical shift"); // NOLINT
  param_manager->declare_parameter(name_prefix +  "LMY"  , &(p_.LMY)   , 1.0 , tam::pmg::ParameterType::DOUBLE, "Scale factor of rolling resistance torque"); // NOLINT
  // clang-format on

  // Other params
  param_manager->declare_parameter(
    name_prefix + "FNOMIN", &(p_.FNOMIN), 3000.0, tam::pmg::ParameterType::DOUBLE,
    "Nominal wheel load");  // NOLINT
}
}  // namespace tam::ocd::tire_models
// MF6.2 steady combined-slip forces and moments, no turn-slip mode (TIR USE_MODE114).
// Equations: Besselink et al., MF6.1 extensions; Pacejka4.E9--4.E75.
// Reference attribution: THIRD_PARTY_LICENSES/MFeval_julia_SOURCE.md.
namespace tam::ocd::tire_models
{
types::TireModelOutput MF52::evaluate62(types::TireModelInput const & u) const
{
  types::TireModelOutput out;
  if (!(std::isfinite(u.inflation_pressure_Pa) && u.inflation_pressure_Pa > 0. && p_.NOMPRES > 0.))
    throw std::invalid_argument("MF62 requires positive gauge inflation pressure");
  if (u.F_z_N <= 0.) return out;
  const double F = u.F_z_N, F0 = p_.LFZO * p_.FNOMIN, d = (F - F0) / F0,
               dp = (u.inflation_pressure_Pa - p_.NOMPRES) / p_.NOMPRES;
  const double a = std::tan(-u.slip_angle_rad), g = u.gamma_rad, gs = std::sin(g),
               k = u.longitudinal_slip;
  auto safe = [](double x) { return x + std::copysign(1e-6, x); };
  auto mf = [](double B, double C, double E, double x) {
    return std::sin(C * std::atan(B * x - E * (B * x - std::atan(B * x))));
  };
  auto weight = [&](double B, double C, double E, double x, double shift) {
    return std::cos(
             C * std::atan(B * (x + shift) - E * (B * (x + shift) - std::atan(B * (x + shift))))) /
           safe(std::cos(C * std::atan(B * shift - E * (B * shift - std::atan(B * shift)))));
  };
  const double Kx = F * (p_.PKX1 + p_.PKX2 * d) * std::exp(p_.PKX3 * d) *
                    (1. + p_.PPX1 * dp + p_.PPX2 * dp * dp) * p_.LKX;
  const double Cx = p_.PCX1 * p_.LCX, Dx = (p_.PDX1 + p_.PDX2 * d) *
                                           (1. + p_.PPX3 * dp + p_.PPX4 * dp * dp) *
                                           (1. - p_.PDX3 * g * g) * p_.LMUX * F;
  const double kx = k + (p_.PHX1 + p_.PHX2 * d) * p_.LHX,
               Ex = std::min(
                 1., (p_.PEX1 + p_.PEX2 * d + p_.PEX3 * d * d) *
                       (1. - p_.PEX4 * std::copysign(1., kx)) * p_.LEX);
  const double fx0 =
    Dx * mf(Kx / safe(Cx * Dx), Cx, Ex, kx) + F * (p_.PVX1 + p_.PVX2 * d) * p_.LVX * p_.LMUX;
  const double Gx = weight(
    (p_.RBX1 + p_.RBX3 * gs * gs) * std::cos(std::atan(p_.RBX2 * k)) * p_.LXAL, p_.RCX1,
    std::min(1., p_.REX1 + p_.REX2 * d), a, p_.RHX1);
  out.longitudinal_force_N = Gx * fx0;
  struct Lat
  {
    double force, K, sh, sv, B, C, mu;
  };
  auto lat = [&](double gg) {
    const double K =
      p_.PKY1 * F0 * (1. + p_.PPY1 * dp) * (1. - p_.PKY3 * std::abs(gg)) *
      std::sin(
        p_.PKY4 * std::atan(F / ((p_.PKY2 + p_.PKY5 * gg * gg) * F0 * (1. + p_.PPY2 * dp)))) *
      p_.LKY;
    const double Kcam = F * (p_.PKY6 + p_.PKY7 * d) * (1. + p_.PPY5 * dp) * p_.LKYC;
    const double svg = F * (p_.PVY3 + p_.PVY4 * d) * gg * p_.LKYC * p_.LMUY;
    const double sv = F * (p_.PVY1 + p_.PVY2 * d) * p_.LVY * p_.LMUY + svg;
    const double sh = (p_.PHY1 + p_.PHY2 * d) * p_.LHY + (Kcam * gg - svg) / safe(K), ay = a + sh,
                 C = p_.PCY1 * p_.LCY;
    const double mu = (p_.PDY1 + p_.PDY2 * d) * (1. + p_.PPY3 * dp + p_.PPY4 * dp * dp) *
                      (1. - p_.PDY3 * gg * gg) * p_.LMUY,
                 D = mu * F, B = K / safe(C * D);
    const double E = std::min(
      1., (p_.PEY1 + p_.PEY2 * d) *
            (1. + p_.PEY5 * gg * gg - (p_.PEY3 + p_.PEY4 * gg) * std::copysign(1., ay)) * p_.LEY);
    return Lat{D * mf(B, C, E, ay) + sv, K, sh, sv, B, C, mu};
  };
  const auto y = lat(gs), y0 = lat(0.);
  const double Gy = weight(
    (p_.RBY1 + p_.RBY4 * gs * gs) * std::cos(std::atan(p_.RBY2 * (a - p_.RBY3))) * p_.LYKA, p_.RCY1,
    std::min(1., p_.REY1 + p_.REY2 * d), k, p_.RHY1 + p_.RHY2 * d);
  const double svyk = y.mu * F * (p_.RVY1 + p_.RVY2 * d + p_.RVY3 * gs) *
                      std::cos(std::atan(p_.RVY4 * a)) *
                      std::sin(p_.RVY5 * std::atan(p_.RVY6 * k)) * p_.LVYKA;
  out.lateral_force_N = Gy * y.force + svyk;
  const double fy = out.lateral_force_N, fx = out.longitudinal_force_N, R = p_.UNLOADED_RADIUS,
               Fm = F < p_.FZMIN ? F * std::pow(F / p_.FZMIN, 2) : F;
  out.overturning_moment_Nm =
    R * Fm * p_.LMX *
      (p_.QSX1 * p_.LVMX - p_.QSX2 * g * (1. + p_.PPMX1 * dp) + p_.QSX3 * fy / p_.FNOMIN +
       p_.QSX4 * std::cos(p_.QSX5 * std::atan(std::pow(p_.QSX6 * Fm / p_.FNOMIN, 2))) *
         std::sin(p_.QSX7 * g + p_.QSX8 * std::atan(p_.QSX9 * fy / p_.FNOMIN)) +
       p_.QSX10 * std::atan(p_.QSX11 * Fm / p_.FNOMIN) * g) +
    R * p_.LMX * (fy * (p_.QSX13 + p_.QSX14 * std::abs(g)) - Fm * p_.QSX12 * g * std::abs(g));
  const double ar = a + y.sh + y.sv / safe(y.K),
               at = a + p_.QHZ1 + p_.QHZ2 * d + (p_.QHZ3 + p_.QHZ4 * d) * gs;
  const double Bt = (p_.QBZ1 + p_.QBZ2 * d + p_.QBZ3 * d * d) *
                    (1. + p_.QBZ4 * g + p_.QBZ5 * std::abs(g)) * p_.LKY / safe(p_.LMUY),
               Ct = p_.QCZ1;
  const double Dt = (p_.QDZ1 + p_.QDZ2 * d) * (1. - p_.PPZ1 * dp) *
                    (1. + p_.QDZ3 * g + p_.QDZ4 * g * g) * F * (R / F0) * p_.LTR;
  const double Et = std::min(
    1., (p_.QEZ1 + p_.QEZ2 * d + p_.QEZ3 * d * d) *
          (1. + (p_.QEZ4 + p_.QEZ5 * gs) * 2. / M_PI * std::atan(Bt * Ct * at)));
  const double Br = p_.QBZ9 * p_.LKY / safe(p_.LMUY) + p_.QBZ10 * y.B * y.C;
  const double Dr =
    F * R *
    ((p_.QDZ6 + p_.QDZ7 * d) * p_.LRES +
     ((p_.QDZ8 + p_.QDZ9 * d) * (1. + p_.PPZ2 * dp) + (p_.QDZ10 + p_.QDZ11 * d) * std::abs(gs)) *
       gs * p_.LKZC) *
    p_.LMUY * std::cos(a);
  const double ratio = std::pow(Kx / safe(y.K) * k, 2);
  const double areq = std::copysign(std::atan(std::sqrt(std::pow(std::tan(ar), 2) + ratio)), ar),
               ateq = std::copysign(std::atan(std::sqrt(std::pow(std::tan(at), 2) + ratio)), at);
  const double trail =
    Dt * std::cos(Ct * std::atan(Bt * ateq - Et * (Bt * ateq - std::atan(Bt * ateq)))) *
    std::cos(u.slip_angle_rad) * p_.LFZO;
  const double arm = R * (p_.SSZ1 + p_.SSZ2 * fy / p_.FNOMIN + (p_.SSZ3 + p_.SSZ4 * d) * g) * p_.LS;
  out.self_aligning_moment_Nm =
    -trail * Gy * y0.force + Dr * std::cos(std::atan(Br * areq)) + arm * fx;
  return out;
}
}  // namespace tam::ocd::tire_models
