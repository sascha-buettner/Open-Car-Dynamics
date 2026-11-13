// Copyright 2026 TUM Autonomous Motorsport
#pragma once
namespace tam::ocd::vehicle_dynamics
{
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calc_dynamic_tire_radius()
{
  // Interpolation function
  imr_.dynamic_tire_radius_m =
    double_per_wheel_t::from_front_and_rear(p_.rr_w_f, p_.rr_w_r) *
    tam::helpers::numerical::interp(
      x_vec_[x::v_x_mps], p_.rr_vel_scale_vel_points_mps, p_.rr_vel_scale_scale_factors);
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<
  TireModelT, AeroModelT>::calculate_tire_spring_initial_compression_m()
{
  if (p_.radial_lookup_enabled > .5) {
    for (int i = 0; i < 4; ++i) {
      const double load = p_.radial_conditioning_enabled > .5 && p_.radial_frozen_preload > .5
                            ? p_.radial_condition_load_N[i]
                            : (p_.corner_weights_enabled > .5
                                 ? p_.corner_weights_kg[i] * G_CONST
                                 : p_.m * G_CONST * (i < 2 ? p_dep_.l_r : p_.l_f) / (2 * p_.l));
      radial_preload_N_[i] = load;
      if (p_.radial_conditioning_enabled > .5)
        imr_.tire_spring_initial_compression_m[i] =
          load / (radial_k_Npm_[i] / p_.radial_stiffness_scale);
      else {
        std::vector<double> signature{
          p_.m, p_.l, p_.l_f, p_.corner_weights_enabled, p_.atmosphere_Pa};
        for (int j = 0; j < 4; ++j) {
          signature.push_back(p_.corner_weights_kg[j]);
          signature.push_back(p_.inflation_reference_absolute_Pa[j]);
        }
        signature.push_back(p_.radial_nominal_stiffness);
        signature.push_back(p_.c_t_f);
        signature.push_back(p_.c_t_r);
        if (signature != radial_compliance_reference_signature_) {
          for (int j = 0; j < 4; ++j) {
            const double f = p_.corner_weights_enabled > .5
                               ? p_.corner_weights_kg[j] * G_CONST
                               : p_.m * G_CONST * (j < 2 ? p_dep_.l_r : p_.l_f) / (2 * p_.l);
            GenericRadialCompliance curve;
            curve.reset(j >= 2, 0., p_.inflation_reference_absolute_Pa[j] - p_.atmosphere_Pa);
            radial_nominal_compression_m_[j] = p_.radial_nominal_stiffness > .5
                                                 ? f / (j < 2 ? p_.c_t_f : p_.c_t_r)
                                                 : curve.values(f)[0];
          }
          radial_compliance_reference_signature_ = signature;
        }
        imr_.tire_spring_initial_compression_m[i] = radial_nominal_compression_m_[i];
      }
    }
    return;
  }
  if (p_.corner_weights_enabled > .5) {
    for (int i = 0; i < 4; ++i)
      imr_.tire_spring_initial_compression_m[i] =
        p_.corner_weights_kg[i] * G_CONST / (i < 2 ? p_.c_t_f : p_.c_t_r);
    return;
  }
  // Tire spring compression @ standstill
  imr_.tire_spring_initial_compression_m.front_left =
    (p_.m * G_CONST * p_dep_.l_r / (2 * p_.l)) / p_.c_t_f;
  imr_.tire_spring_initial_compression_m.rear_left =
    (p_.m * G_CONST * p_.l_f / (2 * p_.l)) / p_.c_t_r;
  imr_.tire_spring_initial_compression_m.front_right =
    imr_.tire_spring_initial_compression_m.front_left;
  imr_.tire_spring_initial_compression_m.rear_right =
    imr_.tire_spring_initial_compression_m.rear_left;
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_tire_spring_force_N()
{  // Tire spring forces
  point_road_gx_ = {0, 0, 0, 0};
  point_road_gy_ = {0, 0, 0, 0};
  point_road_h_ = external_influences_.z_height_road_m;
  if (road_surface_ && p_.contact_enabled < .5) {
    road_surface_->prepare(x_vec_[x::x_m], x_vec_[x::y_m], x_vec_[x::psi_rad]);
    const double c = std::cos(x_vec_[x::psi_rad]), ss = std::sin(x_vec_[x::psi_rad]);
    for (int i = 0; i < 4; ++i) {
      const double wx = i < 2 ? p_.l_f : -p_dep_.l_r,
                   wy = (i % 2 ? -1. : 1.) * (i < 2 ? p_.b_f : p_.b_r) / 2.;
      const auto q =
        road_surface_->sample(x_vec_[x::x_m] + c * wx - ss * wy, x_vec_[x::y_m] + ss * wx + c * wy);
      if (!q.valid) throw std::domain_error("Legacy spring contact leaves diagnostic residual map");
      point_road_h_[i] += q.height_m;
      point_road_gx_[i] = q.gradient_x;
      point_road_gy_[i] = q.gradient_y;
    }
  }

  if (p_.mf62_vertical_enabled > .5) {
    if (p_.belt_enabled > .5)
      throw std::invalid_argument("Use the rolling compliance law for the vertical-belt branch");
    std::vector<double> signature{p_.m, p_.l, p_.l_f};
    for (int i = 0; i < 4; ++i) {
      signature.push_back(p_.inflation_reference_absolute_Pa[i]);
      signature.push_back(p_.camber_out_rad[i]);
    }
    signature.push_back(p_.vertical_front.UNLOADED_RADIUS);
    signature.push_back(p_.vertical_front.WIDTH);
    signature.push_back(p_.vertical_front.ASPECT_RATIO);
    signature.push_back(p_.vertical_front.RIM_RADIUS);
    signature.push_back(p_.vertical_front.LONGVL);
    signature.push_back(p_.vertical_front.FNOMIN);
    signature.push_back(p_.vertical_front.NOMPRES);
    signature.push_back(p_.vertical_front.VERTICAL_STIFFNESS);
    signature.push_back(p_.vertical_front.VERTICAL_DAMPING);
    signature.push_back(p_.vertical_front.Q_RE0);
    signature.push_back(p_.vertical_front.Q_V1);
    signature.push_back(p_.vertical_front.Q_V2);
    signature.push_back(p_.vertical_front.Q_FZ2);
    signature.push_back(p_.vertical_front.Q_FCX);
    signature.push_back(p_.vertical_front.Q_FCY);
    signature.push_back(p_.vertical_front.Q_FCY2);
    signature.push_back(p_.vertical_front.Q_CAM1);
    signature.push_back(p_.vertical_front.Q_CAM2);
    signature.push_back(p_.vertical_front.Q_CAM3);
    signature.push_back(p_.vertical_front.Q_FYS1);
    signature.push_back(p_.vertical_front.Q_FYS2);
    signature.push_back(p_.vertical_front.Q_FYS3);
    signature.push_back(p_.vertical_front.PFZ1);
    signature.push_back(p_.vertical_front.BREFF);
    signature.push_back(p_.vertical_front.DREFF);
    signature.push_back(p_.vertical_front.FREFF);
    signature.push_back(p_.vertical_front.BOTTOM_OFFST);
    signature.push_back(p_.vertical_front.BOTTOM_STIFF);
    signature.push_back(p_.vertical_front.QSY1);
    signature.push_back(p_.vertical_front.QSY2);
    signature.push_back(p_.vertical_front.QSY3);
    signature.push_back(p_.vertical_front.QSY4);
    signature.push_back(p_.vertical_front.QSY5);
    signature.push_back(p_.vertical_front.QSY6);
    signature.push_back(p_.vertical_front.QSY7);
    signature.push_back(p_.vertical_front.QSY8);
    signature.push_back(p_.vertical_rear.UNLOADED_RADIUS);
    signature.push_back(p_.vertical_rear.WIDTH);
    signature.push_back(p_.vertical_rear.ASPECT_RATIO);
    signature.push_back(p_.vertical_rear.RIM_RADIUS);
    signature.push_back(p_.vertical_rear.LONGVL);
    signature.push_back(p_.vertical_rear.FNOMIN);
    signature.push_back(p_.vertical_rear.NOMPRES);
    signature.push_back(p_.vertical_rear.VERTICAL_STIFFNESS);
    signature.push_back(p_.vertical_rear.VERTICAL_DAMPING);
    signature.push_back(p_.vertical_rear.Q_RE0);
    signature.push_back(p_.vertical_rear.Q_V1);
    signature.push_back(p_.vertical_rear.Q_V2);
    signature.push_back(p_.vertical_rear.Q_FZ2);
    signature.push_back(p_.vertical_rear.Q_FCX);
    signature.push_back(p_.vertical_rear.Q_FCY);
    signature.push_back(p_.vertical_rear.Q_FCY2);
    signature.push_back(p_.vertical_rear.Q_CAM1);
    signature.push_back(p_.vertical_rear.Q_CAM2);
    signature.push_back(p_.vertical_rear.Q_CAM3);
    signature.push_back(p_.vertical_rear.Q_FYS1);
    signature.push_back(p_.vertical_rear.Q_FYS2);
    signature.push_back(p_.vertical_rear.Q_FYS3);
    signature.push_back(p_.vertical_rear.PFZ1);
    signature.push_back(p_.vertical_rear.BREFF);
    signature.push_back(p_.vertical_rear.DREFF);
    signature.push_back(p_.vertical_rear.FREFF);
    signature.push_back(p_.vertical_rear.BOTTOM_OFFST);
    signature.push_back(p_.vertical_rear.BOTTOM_STIFF);
    signature.push_back(p_.vertical_rear.QSY1);
    signature.push_back(p_.vertical_rear.QSY2);
    signature.push_back(p_.vertical_rear.QSY3);
    signature.push_back(p_.vertical_rear.QSY4);
    signature.push_back(p_.vertical_rear.QSY5);
    signature.push_back(p_.vertical_rear.QSY6);
    signature.push_back(p_.vertical_rear.QSY7);
    signature.push_back(p_.vertical_rear.QSY8);
    signature.push_back(p_.atmosphere_Pa);
    for (int i = 0; i < 4; ++i) signature.push_back(p_.corner_weights_kg[i]);
    const int zz[4] = {x::z_w_FL_m, x::z_w_FR_m, x::z_w_RL_m, x::z_w_RR_m};
    const int vv[4] = {x::v_z_w_FL_mps, x::v_z_w_FR_mps, x::v_z_w_RL_mps, x::v_z_w_RR_mps};
    if (signature != radial_reference_signature_) {
      for (int i = 0; i < 4; ++i) {
        const auto & t = i < 2 ? p_.vertical_front : p_.vertical_rear;
        const double load = p_.corner_weights_enabled > .5
                              ? p_.corner_weights_kg[i] * G_CONST
                              : p_.m * G_CONST * (i < 2 ? p_.l - p_.l_f : p_.l_f) / (2. * p_.l);
        radial_reference_loaded_m_[i] = mf62_loaded_radius(
          t, load, p_.inflation_reference_absolute_Pa[i] - p_.atmosphere_Pa, 0.,
          -p_.camber_out_rad[i]);
      }
      radial_reference_signature_ = signature;
    }
    // Shared spatial contact uses the same zero-load reference as the
    // manufacturer's nonlinear radial law, rather than Fz divided by a tangent.
    if (p_.contact_enabled > .5) {
      for (int i = 0; i < 4; ++i) {
        const auto & tire = i < 2 ? p_.vertical_front : p_.vertical_rear;
        imr_.tire_spring_initial_compression_m[i] =
          mf62_free_radius(tire, 0.) - radial_reference_loaded_m_[i];
      }
      return;
    }
    for (int i = 0; i < 4; ++i) {
      const auto & t = i < 2 ? p_.vertical_front : p_.vertical_rear;
      const double pressure = p_.inflation_absolute_Pa[i] - p_.atmosphere_Pa;
      const double gamma = -(p_.camber_out_rad[i] + kinematic_camber_rad_[i]);
      tire_gamma_rad_[i] = gamma;
      const double omega = wheel_speeds_radps_[i],
                   loaded = radial_reference_loaded_m_[i] + x_vec_[zz[i]] - point_road_h_[i];
      mf62_free_radius_m_[i] = mf62_free_radius(t, omega);
      mf62_loaded_radius_m_[i] = loaded;
      mf62_deflection_m_[i] = std::max(0., mf62_free_radius_m_[i] - loaded);
      double normal = mf62_vertical_force(t, pressure, omega, gamma, loaded), previous = normal;
      for (int iteration = 0; iteration < 20; ++iteration) {
        double fx = 0., fy = 0.;
        if (std::abs(x_vec_[x::v_x_mps]) > .05 || std::abs(omega) > .1) {
          const int ks[4] = {
            x::kappa_tire_FL, x::kappa_tire_FR, x::kappa_tire_RL, x::kappa_tire_RR};
          const int as[4] = {
            x::alpha_tire_FL, x::alpha_tire_FR, x::alpha_tire_RL, x::alpha_tire_RR};
          const double mirror = i % 2 ? -1. : 1.;
          const auto f = tire_models_[i].evaluate(
            {x_vec_[ks[i]], mirror * x_vec_[as[i]], std::max(normal, .001), gamma, pressure});
          fx = f.longitudinal_force_N * external_influences_.lambda_mue[i];
          fy = f.lateral_force_N * external_influences_.lambda_mue[i];
        }
        previous = normal;
        const double next = mf62_vertical_force(t, pressure, omega, gamma, loaded, fx, fy);
        normal = .5 * (normal + next);
        if (std::abs(next - previous) < 1e-5 * std::max(1., normal)) break;
      }
      radial_residual_N_[i] = normal - previous;
      const double wx = i < 2 ? p_.l_f : -p_dep_.l_r,
                   wy = (i % 2 ? -1. : 1.) * (i < 2 ? p_.b_f : p_.b_r) / 2.;
      const double vx = x_vec_[x::v_x_mps] - x_vec_[x::psi_dot_radps] * wy,
                   vy = x_vec_[x::v_y_mps] + x_vec_[x::psi_dot_radps] * wx,
                   c = std::cos(x_vec_[x::psi_rad]), s = std::sin(x_vec_[x::psi_rad]);
      const double road_v =
        point_road_gx_[i] * (c * vx - s * vy) + point_road_gy_[i] * (s * vx + c * vy);
      const double damping = normal > 0. ? -t.VERTICAL_DAMPING * (x_vec_[vv[i]] - road_v) : 0.;
      imr_.tire_spring_force_N[i] = -std::max(0., normal + damping);
    }
    return;
  }

  // Clamp to Fz = 0.001 in order to avoid numeric problems in the tire model.
  imr_.tire_spring_force_N.front_left = std::clamp(
    (p_.radial_lookup_enabled > .5 ? radial_force_coefficient_Npm_[0] : p_.c_t_f) *
      (x_vec_[x::z_w_FL_m] - point_road_h_.front_left -
       imr_.tire_spring_initial_compression_m.front_left),
    std::numeric_limits<double>::lowest(), -0.001);
  imr_.tire_spring_force_N.front_right = std::clamp(
    (p_.radial_lookup_enabled > .5 ? radial_force_coefficient_Npm_[1] : p_.c_t_f) *
      (x_vec_[x::z_w_FR_m] - point_road_h_.front_right -
       imr_.tire_spring_initial_compression_m.front_right),
    std::numeric_limits<double>::lowest(), -0.001);
  imr_.tire_spring_force_N.rear_left = std::clamp(
    (p_.radial_lookup_enabled > .5 ? radial_force_coefficient_Npm_[2] : p_.c_t_r) *
      (x_vec_[x::z_w_RL_m] - point_road_h_.rear_left -
       imr_.tire_spring_initial_compression_m.rear_left),
    std::numeric_limits<double>::lowest(), -0.001);
  imr_.tire_spring_force_N.rear_right = std::clamp(
    (p_.radial_lookup_enabled > .5 ? radial_force_coefficient_Npm_[3] : p_.c_t_r) *
      (x_vec_[x::z_w_RR_m] - point_road_h_.rear_right -
       imr_.tire_spring_initial_compression_m.rear_right),
    std::numeric_limits<double>::lowest(), -0.001);
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_local_road_contact()
{
  imr_.contact_force_x_N = {0, 0, 0, 0};
  imr_.contact_force_y_N = {0, 0, 0, 0};
  imr_.contact_missing_samples = {0, 0, 0, 0};
  imr_.contact_bridged_samples = {0, 0, 0, 0};
  imr_.contact_yaw_moment_Nm = {0, 0, 0, 0};
  imr_.contact_elastic_energy_J = {0, 0, 0, 0};
  imr_.tire_deformation_extrapolated = {0, 0, 0, 0};
  imr_.road_contact_height_m = external_influences_.z_height_road_m;
  for (int node = 0; node < 5; ++node) {
    imr_.contact_node_query_x_m[node] = {0, 0, 0, 0};
    imr_.contact_node_query_y_m[node] = {0, 0, 0, 0};
    imr_.contact_node_road_z_m[node] = {0, 0, 0, 0};
    imr_.contact_node_belt_z_m[node] = {0, 0, 0, 0};
    imr_.contact_node_normal_N[node] = {0, 0, 0, 0};
    imr_.contact_node_force_x_N[node] = {0, 0, 0, 0};
    imr_.contact_node_force_y_N[node] = {0, 0, 0, 0};
  }
  imr_.tire_lateral_drop_m = {0, 0, 0, 0};
  if (p_.belt_enabled > .5 && p_.contact_enabled < .5)
    throw std::invalid_argument("Dynamic tire belt requires five-point contact");
  if (p_.contact_enabled < 0.5) {
    if (road_surface_)
      for (int i = 0; i < 4; ++i) {
        imr_.road_contact_height_m[i] = point_road_h_[i];
        const double normal = -imr_.tire_spring_force_N[i], c = std::cos(x_vec_[x::psi_rad]),
                     ss = std::sin(x_vec_[x::psi_rad]);
        imr_.contact_force_x_N[i] = -normal * (c * point_road_gx_[i] + ss * point_road_gy_[i]);
        imr_.contact_force_y_N[i] = -normal * (-ss * point_road_gx_[i] + c * point_road_gy_[i]);
        imr_.contact_elastic_energy_J[i] =
          .5 * normal * normal /
          (p_.radial_lookup_enabled > .5 ? radial_force_coefficient_Npm_[i]
                                         : (i < 2 ? p_.c_t_f : p_.c_t_r));
      }
    return;
  }
  if (p_.contact_sampling_points != 1. && p_.contact_sampling_points != 5.)
    throw std::invalid_argument("Contact sampling_points must be exactly 1 or 5");
  for (double positive :
       {p_.c_t_f, p_.c_t_r, p_.lateral_stiffness_f, p_.lateral_stiffness_r, p_.sidewall_height_f,
        p_.sidewall_height_r, p_.tire_width_f, p_.tire_width_r})
    if (!(std::isfinite(positive) && positive > 0.))
      throw std::invalid_argument(
        "Contact stiffness, sidewall height and tire width must be positive and finite");
  for (double nonnegative :
       {p_.contact_half_length, p_.contact_micro_radius, p_.enveloping_radius_f,
        p_.enveloping_radius_r, p_.tire_damping_f, p_.tire_damping_r, p_.lateral_drop_scale})
    if (!(std::isfinite(nonnegative) && nonnegative >= 0.))
      throw std::invalid_argument(
        "Contact geometry, damping and lateral coupling must be nonnegative and finite");
  if (road_surface_) road_surface_->prepare(x_vec_[x::x_m], x_vec_[x::y_m], x_vec_[x::psi_rad]);
  const int latidx[4] = {
    x::lateral_deflection_tire_FL_m, x::lateral_deflection_tire_FR_m,
    x::lateral_deflection_tire_RL_m, x::lateral_deflection_tire_RR_m};
  const int zidx[4] = {
    p_.belt_enabled > .5 ? x::z_belt_FL_m : x::z_w_FL_m,
    p_.belt_enabled > .5 ? x::z_belt_FR_m : x::z_w_FR_m,
    p_.belt_enabled > .5 ? x::z_belt_RL_m : x::z_w_RL_m,
    p_.belt_enabled > .5 ? x::z_belt_RR_m : x::z_w_RR_m};
  const int vidx[4] = {
    p_.belt_enabled > .5 ? x::v_z_belt_FL_mps : x::v_z_w_FL_mps,
    p_.belt_enabled > .5 ? x::v_z_belt_FR_mps : x::v_z_w_FR_mps,
    p_.belt_enabled > .5 ? x::v_z_belt_RL_mps : x::v_z_w_RL_mps,
    p_.belt_enabled > .5 ? x::v_z_belt_RR_mps : x::v_z_w_RR_mps};
  const double wx[4] = {p_.l_f, p_.l_f, -p_dep_.l_r, -p_dep_.l_r};
  const double wy[4] = {p_.b_f / 2, -p_.b_f / 2, p_.b_r / 2, -p_.b_r / 2};
  const double cp = std::cos(x_vec_[x::psi_rad]), sp = std::sin(x_vec_[x::psi_rad]);
  for (int w = 0; w < 4; ++w) {
    const double delta = imr_.effective_steering_angle_per_wheel_rad[w];
    const double cd = std::cos(delta), sd = std::sin(delta);
    const double total_kt = p_.radial_lookup_enabled > .5 ? radial_force_coefficient_Npm_[w]
                                                          : (w < 2 ? p_.c_t_f : p_.c_t_r);
    double kt = total_kt, damping = w < 2 ? p_.tire_damping_f : p_.tire_damping_r;
    if (p_.belt_enabled > .5) {
      const double mass = w < 2 ? p_.belt_mass_f : p_.belt_mass_r;
      const double frequency = w < 2 ? p_.belt_frequency_f : p_.belt_frequency_r;
      const double ratio = w < 2 ? p_.belt_damping_ratio_f : p_.belt_damping_ratio_r;
      const double cb = mass * std::pow(2. * M_PI * frequency, 2);
      if (!(std::isfinite(mass) && std::isfinite(frequency) && std::isfinite(ratio) && mass > 0. &&
            mass < (w < 2 ? p_.m_t_f : p_.m_t_r) && frequency > 0. && ratio >= 0. && cb > total_kt))
        throw std::invalid_argument(
          "Belt mass/frequency must leave positive rim mass and residual compliance");
      kt = 1. / (1. / total_kt - 1. / cb);
      const double db = 2. * ratio * mass * 2. * M_PI * frequency;
      damping =
        std::max(0., (damping - db * std::pow(total_kt / cb, 2)) / std::pow(total_kt / kt, 2));
      if (p_.radial_lookup_enabled > .5 && p_.radial_conditioning_enabled < .5) {
        const double tangent = radial_k_Npm_[w], contact_tangent = 1. / (1. / tangent - 1. / cb);
        damping = std::max(
          0., ((w < 2 ? p_.tire_damping_f : p_.tire_damping_r) - db * std::pow(tangent / cb, 2)) /
                std::pow(tangent / contact_tangent, 2));
      }
    }
    const double width = w < 2 ? p_.tire_width_f : p_.tire_width_r;
    const double lateral = p_.sidewall_enabled > .5 ? x_vec_[latidx[w]] : 0.;
    const double sidewall = w < 2 ? p_.sidewall_height_f : p_.sidewall_height_r;
    if (!(std::isfinite(p_.lateral_drop_deflection_limit) &&
          p_.lateral_drop_deflection_limit >= 0. && std::isfinite(sidewall) && sidewall > 0.))
      throw std::invalid_argument("Lateral height-coupling bound/sidewall height invalid");
    const double limited =
      std::clamp(lateral, -p_.lateral_drop_deflection_limit, p_.lateral_drop_deflection_limit);
    const double ratio = std::clamp(limited / sidewall, -.95, .95);
    const double drop = p_.lateral_drop_scale * sidewall * (1 - std::sqrt(1 - ratio * ratio));
    imr_.tire_lateral_deflection_m[w] = lateral;
    imr_.tire_lateral_drop_m[w] = drop;
    imr_.tire_deformation_extrapolated[w] = std::abs(lateral) > .02 ? 1. : 0.;
    const double vx = x_vec_[x::v_x_mps] - wy[w] * x_vec_[x::psi_dot_radps];
    const double vy = x_vec_[x::v_y_mps] + wx[w] * x_vec_[x::psi_dot_radps];
    const double vxglobal = cp * vx - sp * vy, vyglobal = sp * vx + cp * vy;
    std::array<double, 5> heights{}, rates{}, gx{}, gy{}, heading_derivative{};
    // Collapsing the massless samples to their centre supplies a point contact
    // with exactly the same constitutive law, belt and transient states.
    const double length = p_.contact_sampling_points == 5. ? p_.contact_half_length : 0.;
    const double half_width = p_.contact_sampling_points == 5. ? width / 2 : 0.;
    const double px[5] = {0., length, length, -length, -length};
    const double py[5] = {
      lateral, half_width + lateral, -half_width + lateral, half_width + lateral,
      -half_width + lateral};
    for (int i = 0; i < 5; ++i) {
      const double bx = wx[w] + cd * px[i] - sd * py[i], by = wy[w] + sd * px[i] + cd * py[i];
      const double qx = x_vec_[x::x_m] + cp * bx - sp * by, qy = x_vec_[x::y_m] + sp * bx + cp * by;
      double h = external_influences_.z_height_road_m[w], dx = 0, dy = 0;
      if (road_surface_) {
        const auto sample = road_surface_->tire_contact(
          qx, qy, x_vec_[x::psi_rad] + delta,
          w < 2 ? p_.enveloping_radius_f : p_.enveloping_radius_r, p_.contact_micro_radius);
        if (sample.valid) {
          h += sample.height_m;
          if (sample.bridged_missing_centre) imr_.contact_bridged_samples[w] += 1.;
          dx = sample.gradient_x;
          dy = sample.gradient_y;
          heading_derivative[i] = sample.heading_derivative_m;
        } else {
          imr_.contact_missing_samples[w] += 1;
          if (p_.strict_surface_coverage > .5)
            throw std::domain_error(
              "Tire collider leaves measured surface coverage: wheel=" + std::to_string(w) +
              " node=" + std::to_string(i) + " planar_x_m=" + std::to_string(qx) + " planar_y_m=" +
              std::to_string(qy) + "; inspect contact_missing_samples and the map crop");
        }
      }
      imr_.contact_node_query_x_m[i][w] = qx;
      imr_.contact_node_query_y_m[i][w] = qy;
      imr_.contact_node_road_z_m[i][w] = h;
      heights[i] = h + imr_.tire_spring_initial_compression_m[w] * (total_kt / kt) -
                   x_vec_[zidx[w]] - (p_.belt_enabled > .5 ? 0. : drop);
      if (p_.radial_lookup_enabled > .5 && p_.radial_conditioning_enabled < .5) {
        double reference = imr_.tire_spring_initial_compression_m[w];
        if (p_.belt_enabled > .5) {
          const double mass = w < 2 ? p_.belt_mass_f : p_.belt_mass_r,
                       frequency = w < 2 ? p_.belt_frequency_f : p_.belt_frequency_r;
          reference -= radial_preload_N_[w] / (mass * std::pow(2. * M_PI * frequency, 2));
        }
        heights[i] = h + reference - x_vec_[zidx[w]] - (p_.belt_enabled > .5 ? 0. : drop);
      }
      const double pointvx =
        vxglobal + x_vec_[x::psi_dot_radps] * (-sp * (bx - wx[w]) - cp * (by - wy[w]));
      const double pointvy =
        vyglobal + x_vec_[x::psi_dot_radps] * (cp * (bx - wx[w]) - sp * (by - wy[w]));
      rates[i] = dx * pointvx + dy * pointvy - x_vec_[vidx[w]] +
                 heading_derivative[i] * x_vec_[x::psi_dot_radps];
      gx[i] = cp * dx + sp * dy;
      gy[i] = -sp * dx + cp * dy;
      if (i == 0) imr_.road_contact_height_m[w] = h;
    }
    auto contact = contact_.evaluate(heights, rates, kt, damping);
    if (p_.mf62_vertical_enabled > .5) {
      const auto & tire = w < 2 ? p_.vertical_front : p_.vertical_rear;
      const double omega = wheel_speeds_radps_[w];
      const double growth = mf62_free_radius(tire, omega) - mf62_free_radius(tire, 0.);
      auto nonlinear_heights = heights;
      for (double & value : nonlinear_heights) value += growth;
      const int ks[4] = {x::kappa_tire_FL, x::kappa_tire_FR, x::kappa_tire_RL, x::kappa_tire_RR};
      const int as[4] = {x::alpha_tire_FL, x::alpha_tire_FR, x::alpha_tire_RL, x::alpha_tire_RR};
      const double pressure = p_.inflation_absolute_Pa[w] - p_.atmosphere_Pa;
      const double gamma = -(p_.camber_out_rad[w] + kinematic_camber_rad_[w]);
      auto tangential = [&](double normal) {
        if (normal <= 0.) return std::array<double, 2>{0., 0.};
        const double load = std::max(normal, p_.low_load_boundary_N);
        const double mirror = w % 2 ? -1. : 1.;
        const double kappa =
          p_.force_conditioning_enabled > .5 ? p_.force_condition_kappa[w] : x_vec_[ks[w]];
        const double alpha =
          p_.force_conditioning_enabled > .5 ? p_.force_condition_alpha[w] : x_vec_[as[w]];
        const auto force = tire_models_[w].evaluate({kappa, mirror * alpha, load, gamma, pressure});
        const double memory =
          p_.load_memory_fraction > 0.
            ? load_capacity_factor(
                normal, x_vec_[x::contact_capacity_FL_N + w], p_.load_memory_fraction)
            : 1.;
        const double scale = normal / load * memory * external_influences_.lambda_mue[w];
        return std::array<double, 2>{
          scale * force.longitudinal_force_N, scale * force.lateral_force_N};
      };
      contact = nonlinear_radial_contact(
        contact_, nonlinear_heights, rates, tire, pressure, omega, gamma,
        p_.mf62_contact_force_coupling > .5, tangential);
      mf62_free_radius_m_[w] = mf62_free_radius(tire, omega);
      mf62_deflection_m_[w] =
        std::sqrt(2. * contact_.evaluate(nonlinear_heights, rates, 1., 0.).elastic_energy_J);
      mf62_loaded_radius_m_[w] = mf62_free_radius_m_[w] - mf62_deflection_m_[w];
    }
    imr_.tire_spring_force_N[w] = -contact.force_N;
    imr_.contact_elastic_energy_J[w] = contact.elastic_energy_J;
    for (int i = 0; i < 5; ++i) {
      imr_.contact_node_normal_N[i][w] = contact.normal_N[i];
      imr_.contact_node_force_x_N[i][w] = -contact.normal_N[i] * gx[i];
      imr_.contact_node_force_y_N[i][w] = -contact.normal_N[i] * gy[i];
      // h - compression + elastic displacement: active nodes lie on the
      // enveloped support surface; free nodes can lift clear of it.
      imr_.contact_node_belt_z_m[i][w] =
        imr_.contact_node_road_z_m[i][w] - heights[i] + contact.belt_displacement_m[i];
      imr_.contact_force_x_N[w] -= contact.normal_N[i] * gx[i];
      imr_.contact_force_y_N[w] -= contact.normal_N[i] * gy[i];
      const double bx = cd * px[i] - sd * py[i], by = sd * px[i] + cd * py[i];
      imr_.contact_yaw_moment_Nm[w] +=
        contact.normal_N[i] * (by * gx[i] - bx * gy[i] - heading_derivative[i]);
    }
  }
}
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_vertical_tire_force_N()
{
  for (int i = 0; i < 4; ++i) {
    const double vertical = -imr_.tire_spring_force_N[i];
    imr_.vertical_tire_force_N[i] = vertical;
    double gx = 0., gy = 0.;
    if (vertical > 1e-8) {
      const double delta = imr_.effective_steering_angle_per_wheel_rad[i];
      const double body_x = -imr_.contact_force_x_N[i] / vertical;
      const double body_y = -imr_.contact_force_y_N[i] / vertical;
      gx = std::cos(delta) * body_x + std::sin(delta) * body_y;
      gy = -std::sin(delta) * body_x + std::cos(delta) * body_y;
    }
    contact_frames_[i] = RoadContactFrame(gx, gy);
    imr_.contact_normal_load_N[i] = vertical * contact_frames_[i].normal_scale;
  }
}
}  // namespace tam::ocd::vehicle_dynamics
