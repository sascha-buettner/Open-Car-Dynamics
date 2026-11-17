// Copyright 2026 Sascha Buettner
#pragma once
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "ocd_point_mass_model_cpp/dual2.hpp"
#include "ocd_point_mass_model_cpp/smooth_table.hpp"

namespace tam::ocd::point_mass
{
constexpr double kGravity = 9.81;

/// Track quantities at one arclength (constant per collocation node): slope mu, banking phi and
/// the angular rates Omega of the reference frame per metre with their derivatives in s.
struct TrackPoint
{
  double mu{0.0}, phi{0.0};
  double Omega_x{0.0}, dOmega_x{0.0};
  double Omega_y{0.0}, dOmega_y{0.0};
  double Omega_z{0.0}, dOmega_z{0.0};
};

/// Parameters of the point mass on the 3D track (the raceline optimisation's cgr keys).
struct PointMassParams
{
  double h_cog_m{0.0};  // vehicle_params.h
  double Iz_kgm2{0.0};  // vehicle_params.Iz
  // objective: w_T / s_dot + jerk, curvature and late-apex regularisation
  double w_T{1.0}, w_jerk_x{0.0}, w_jerk_y{0.0}, w_jerk_y_gain{0.0}, w_jerk_xy{0.0};
  double w_kappa{0.0}, w_late_apex{0.0};
  bool curvature_penalty_sdot_scaled{true};
  // simplifications of the apparent accelerations
  bool neglect_w_omega_y{false}, neglect_w_omega_x{false}, neglect_euler{false};
  bool neglect_centrifugal{false}, neglect_w_dot{false}, neglect_V_omega{false};
  // diamond limit, power_rows formulation
  double abs_epsilon_mps2{1.0e-3};  // smooth |ax|, |ay| inside the powers
  double gg_exponent_ax_pos{2.0}, gg_exponent_ax_neg{2.0};
  bool yaw_moment_limit{false};  // rows on the yaw moment (vehicles with yaw-moment data)
};

/// Smooth GGGV limits of the diamond: (V, g_tilde) -> ax_min, ax_max, ay_max, ym_max (scales of
/// the planner applied), and the machine limit ax_machine(V).
struct DiamondLimits
{
  enum Table : int { kAxMin = 0, kAxMax = 1, kAyMax = 2, kYmMax = 3 };
  SmoothTables2D acc;
  SmoothTable1D machine;
};

/// Apparent accelerations of the point mass in the road frame (planning_common track3D
/// calc_apparent_accelerations).
template <class T>
struct Apparent
{
  T ax, ay, g, omega_z_dot;
};

/// Point mass on the 3D track in the curvilinear frame, s as the independent variable
/// (mod_planning raceline_opt, gg_mode diamond, tire_utilization_model power_rows):
///
///   x = (V, n, chi, ax, ay),  u = (jx, jy)
///   dynamics      dx/ds and the integrands q_t (time) and q_reg (regularisation) of the objective
///   path rows     traction, braking (each +- yaw moment), yaw (+-), machine; each >= 0
///
/// All functions are templated on the scalar type: double for values, Dual2<7> for the exact
/// Jacobian and Hessian in (x, u).
class PointMassDiamond
{
public:
  static constexpr int nx = 5;
  static constexpr int nu = 2;
  static constexpr int nz = nx + nu;

  PointMassDiamond(PointMassParams params, DiamondLimits limits)
  : p_(std::move(params)), lim_(std::move(limits))
  {
    if (p_.yaw_moment_limit && lim_.acc.count() < 4)
      throw std::invalid_argument("PointMassDiamond: the yaw-moment limit needs the ym_max table");
  }

  const PointMassParams & params() const { return p_; }
  const DiamondLimits & limits() const { return lim_; }
  /// Number of path rows.
  int n_path() const { return p_.yaw_moment_limit ? 7 : 3; }

  /// dx/ds, q_t and q_reg.
  template <class T>
  void dynamics(const T * x, const T * u, const TrackPoint & t, T * dx, T & q_t, T & q_reg) const
  {
    using std::cos;
    using std::sin;
    const T & V = x[0];
    const T & n = x[1];
    const T & chi = x[2];
    const T & ax = x[3];
    const T & ay = x[4];
    const T & jx = u[0];
    const T & jy = u[1];
    const T cos_chi = cos(chi), sin_chi = sin(chi);
    const T s_dot = V * cos_chi / (1.0 - n * t.Omega_z);
    const T inv_s_dot = 1.0 / s_dot;
    const T w = n * t.Omega_x * s_dot;
    T dV = inv_s_dot * ax;
    if (!p_.neglect_w_omega_y) dV += w * (t.Omega_x * sin_chi - t.Omega_y * cos_chi);
    T dchi = inv_s_dot * ay / V - t.Omega_z;
    if (!p_.neglect_w_omega_x) dchi += w * (t.Omega_x * cos_chi + t.Omega_y * sin_chi) / V;
    dx[0] = dV;
    dx[1] = inv_s_dot * V * sin_chi;
    dx[2] = dchi;
    dx[3] = jx / s_dot;
    dx[4] = jy / s_dot;

    q_t = p_.w_T * inv_s_dot;
    const T jx_s = jx / s_dot, jy_s = jy / s_dot;
    const T Omega_z_rl = t.Omega_z + dchi;
    T reg(0.0);
    if (p_.w_jerk_x != 0.0) reg += p_.w_jerk_x * square(jx_s);
    if (p_.w_jerk_y != 0.0) reg += p_.w_jerk_y * square(jy_s);
    if (p_.w_kappa != 0.0) {
      if (p_.curvature_penalty_sdot_scaled)
        reg += p_.w_kappa * square(Omega_z_rl) * square(s_dot);
      else
        reg += p_.w_kappa * square(Omega_z_rl);
    }
    if (p_.w_jerk_y_gain != 0.0) reg += p_.w_jerk_y_gain * square(jy_s) * square(ay);
    if (p_.w_jerk_xy != 0.0) reg += p_.w_jerk_xy * square(jx_s) * square(jy_s);
    if (p_.w_late_apex != 0.0) reg += p_.w_late_apex * ax * square(Omega_z_rl) * inv_s_dot;
    q_reg = reg;
  }

  /// Apparent accelerations ax~, ay~, g~ and the yaw acceleration omega_z_dot~.
  template <class T>
  Apparent<T> apparent(const T * x, const T * u, const TrackPoint & t) const
  {
    using std::cos;
    using std::fmax;
    using std::sin;
    const T & V = x[0];
    const T & n = x[1];
    const T & chi = x[2];
    const T & ax = x[3];
    const T & ay = x[4];
    const T & jy = u[1];
    const double Ox = t.Omega_x, Oy = t.Omega_y, Oz = t.Omega_z;
    const double dOx = t.dOmega_x, dOy = t.dOmega_y, dOz = t.dOmega_z;
    const T cos_chi = cos(chi), sin_chi = sin(chi);

    const T s_dot = V * cos_chi / (1.0 - n * Oz);
    const T w = n * Ox * s_dot;
    T V_dot = ax;
    if (!p_.neglect_w_omega_y) V_dot += w * (Ox * sin_chi - Oy * cos_chi) * s_dot;
    const T n_dot = V * sin_chi;
    T chi_dot = ay / V - Oz * s_dot;
    if (!p_.neglect_w_omega_x) chi_dot += w * (Ox * cos_chi + Oy * sin_chi) * s_dot / V;
    const T s_ddot =
      ((V_dot * cos_chi - V * sin_chi * chi_dot) * (1.0 - n * Oz) -
       (V * cos_chi) * (-n_dot * Oz - n * dOz * s_dot)) /
      (1.0 + 2.0 * n * Oz + square(n) * (Oz * Oz));

    T omega_x_dot(0.0), omega_y_dot(0.0);
    if (!p_.neglect_euler) {
      omega_x_dot = (dOx * s_dot * cos_chi - Ox * sin_chi * chi_dot + dOy * s_dot * sin_chi +
                     Oy * cos_chi * chi_dot) *
                      s_dot +
                    (Ox * cos_chi + Oy * sin_chi) * s_ddot;
      omega_y_dot = (-dOx * s_dot * sin_chi - Ox * cos_chi * chi_dot + dOy * s_dot * cos_chi -
                     Oy * sin_chi * chi_dot) *
                      s_dot +
                    (-Ox * sin_chi + Oy * cos_chi) * s_ddot;
    }
    T omega_x(0.0), omega_y(0.0), omega_z(0.0);
    if (!p_.neglect_centrifugal) {
      omega_x = (Ox * cos_chi + Oy * sin_chi) * s_dot;
      omega_y = (-Ox * sin_chi + Oy * cos_chi) * s_dot;
      omega_z = Oz * s_dot + chi_dot;
    }
    T w_dot(0.0);
    if (!p_.neglect_w_dot) w_dot = n_dot * Ox * s_dot + n * dOx * square(s_dot) + n * Ox * s_ddot;
    T V_omega(0.0);
    if (!p_.neglect_V_omega) V_omega = (-Ox * sin_chi + Oy * cos_chi) * s_dot * V;

    // yaw acceleration of the vehicle (yaw-moment limit)
    const T inv_V = 1.0 / V;
    const T chi_ddot = jy * inv_V - ay * square(inv_V) * V_dot - dOz * square(s_dot) -
                       Oz * s_ddot +
                       w_dot * (s_dot * inv_V) * (cos_chi * Ox + sin_chi * Oy) +
                       w * (s_dot * inv_V) *
                         (-sin_chi * chi_dot * Ox + cos_chi * dOx * s_dot +
                          cos_chi * chi_dot * Oy + sin_chi * dOy * s_dot) +
                       w * (s_ddot * inv_V) * (cos_chi * Ox + sin_chi * Oy) +
                       w * (s_dot * square(inv_V)) * V_dot * (-cos_chi * Ox - sin_chi * Oy);
    Apparent<T> a;
    a.omega_z_dot = dOz * square(s_dot) + Oz * s_ddot + chi_ddot;

    const double h = p_.h_cog_m;
    const double cos_mu = std::cos(t.mu), sin_mu = std::sin(t.mu), sin_phi = std::sin(t.phi);
    const double cos_phi = std::cos(t.phi);
    a.ax = ax + omega_y_dot * h - omega_z * omega_x * h +
           kGravity * (-sin_mu * cos_chi + cos_mu * sin_phi * sin_chi);
    a.ay = ay + omega_x_dot * h + omega_z * omega_y * h +
           kGravity * (sin_mu * sin_chi + cos_mu * sin_phi * cos_chi);
    a.g = fmax(
      w_dot - V_omega + (square(omega_x) - square(omega_y)) * h + kGravity * cos_mu * cos_phi, 0.0);
    return a;
  }

  /// Path rows (each >= 0) in the order traction (+Mz, -Mz), braking (+Mz, -Mz), yaw (+, -),
  /// machine; without the yaw-moment limit traction, braking, machine.
  template <class T>
  void path(const T * x, const T * u, const TrackPoint & t, T * rows) const
  {
    using std::fabs;
    using std::pow;
    using std::sqrt;
    const Apparent<T> a = apparent(x, u, t);
    const T & V = x[0];
    T ax_min, ax_max, ay_max, ym_max, ax_machine;
    lookup(V, a.g, ax_min, ax_max, ay_max, ym_max, ax_machine);

    const double eps = std::max(p_.abs_epsilon_mps2, 1.0e-12);
    const T root = sqrt(square(a.ax) + eps * eps);
    const T ay_abs = sqrt(square(a.ay) + eps * eps);
    const T Mz = a.omega_z_dot * p_.Iz_kgm2;
    int r = 0;
    const T ax_parts[2] = {0.5 * (root + a.ax), 0.5 * (root - a.ax)};
    const T ax_lims[2] = {ax_max, fabs(ax_min)};
    const double exps[2] = {p_.gg_exponent_ax_pos, p_.gg_exponent_ax_neg};
    for (int side = 0; side < 2; ++side) {
      const double pe = exps[side];
      const T use =
        pow(pow(ax_parts[side] / ax_lims[side], pe) + pow(ay_abs / ay_max, pe), 1.0 / pe);
      if (p_.yaw_moment_limit) {
        rows[r++] = 1.0 - use - Mz / ym_max;
        rows[r++] = 1.0 - use + Mz / ym_max;
      } else {
        rows[r++] = 1.0 - use;
      }
    }
    if (p_.yaw_moment_limit) {
      rows[r++] = 0.5 - Mz / ym_max;
      rows[r++] = 0.5 + Mz / ym_max;
    }
    rows[r++] = ax_machine - a.ax;
  }

private:
  // limit lookups at (V, g): values for double, exact second-order chain rule for Dual2
  void lookup(
    double V, double g, double & ax_min, double & ax_max, double & ay_max, double & ym_max,
    double & ax_machine) const
  {
    if (p_.yaw_moment_limit) {
      std::array<double, 4> out;
      lim_.acc.values(V, g, std::array<int, 4>{0, 1, 2, 3}, out);
      ax_min = out[0], ax_max = out[1], ay_max = out[2], ym_max = out[3];
    } else {
      std::array<double, 3> out;
      lim_.acc.values(V, g, std::array<int, 3>{0, 1, 2}, out);
      ax_min = out[0], ax_max = out[1], ay_max = out[2], ym_max = 0.0;
    }
    ax_machine = lim_.machine(V);
  }

  template <int N>
  void lookup(
    const Dual2<N> & V, const Dual2<N> & g, Dual2<N> & ax_min, Dual2<N> & ax_max,
    Dual2<N> & ay_max, Dual2<N> & ym_max, Dual2<N> & ax_machine) const
  {
    auto apply = [&](const SmoothTables2D::Derivatives & d) {
      return chain2(V, g, d.f, d.fx, d.fy, d.fxx, d.fxy, d.fyy);
    };
    if (p_.yaw_moment_limit) {
      std::array<SmoothTables2D::Derivatives, 4> out;
      lim_.acc.derivatives(V.v, g.v, std::array<int, 4>{0, 1, 2, 3}, out);
      ax_min = apply(out[0]), ax_max = apply(out[1]), ay_max = apply(out[2]);
      ym_max = apply(out[3]);
    } else {
      std::array<SmoothTables2D::Derivatives, 3> out;
      lim_.acc.derivatives(V.v, g.v, std::array<int, 3>{0, 1, 2}, out);
      ax_min = apply(out[0]), ax_max = apply(out[1]), ay_max = apply(out[2]);
      ym_max = Dual2<N>(0.0);
    }
    double f, f1, f2;
    lim_.machine.eval(V.v, f, f1, f2);
    ax_machine = chain(V, f, f1, f2);
  }

  PointMassParams p_;
  DiamondLimits lim_;
};

}  // namespace tam::ocd::point_mass
