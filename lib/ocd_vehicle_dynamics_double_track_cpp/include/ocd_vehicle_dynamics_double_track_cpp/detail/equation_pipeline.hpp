// Copyright 2026 TUM Autonomous Motorsport
#pragma once
namespace tam::ocd::vehicle_dynamics
{
template <
  tam::ocd::interfaces::concepts::TireModel TireModelT,
  tam::ocd::interfaces::concepts::AerodynamicsModel AeroModelT>
void VehicleDynamicsDoubleTrackEqns<TireModelT, AeroModelT>::calculate_intermediate_results()
{
  calculate_dependent_parameters();
  calculate_effective_steering_angle();
  if(p_.ring_enabled>.5)for(int i=0;i<4;++i){
    imr_.effective_steering_angle_per_wheel_rad[i]+=x_vec_[x::ring_yaw_FL_rad+i];
    kinematic_camber_rad_[i]=(p_.kinematics_enabled>.5?kinematic_camber_rad_[i]:0.)+(i%2?1.:-1.)*x_vec_[x::ring_camber_FL_rad+i];
  }
  calc_dynamic_tire_radius();
  // Iterate only the radial constitutive law with the current ODE-stage forces.
  // No previous-stage forces, prescribed tangential forces, or new states are used.
  if (p_.radial_lookup_enabled > .5) {
    if (!radial_lookup_.configured())
      radial_lookup_.configure(
        p_.radial_grid_speeds, p_.radial_grid_pressures, p_.radial_grid_loads,
        p_.radial_grid_front, p_.radial_grid_rear);
    if (!(std::isfinite(p_.radial_stiffness_scale) && p_.radial_stiffness_scale > 0.))
      throw std::invalid_argument("Invalid radial stiffness scale");
    if (p_.mf62_vertical_enabled > .5)
      throw std::invalid_argument("Select lookup or MF62 radial law, not both");
    for (int i = 0; i < 4; ++i) {
      const double v = radial_lookup_speed_mps();
      const double f = p_.radial_conditioning_enabled > .5 ? p_.radial_condition_load_N[i]
                                                           : p_.corner_weights_kg[i] * G_CONST;
      const double fx = p_.radial_conditioning_enabled > .5 ? p_.radial_condition_fx_N[i] : 0.;
      radial_k_Npm_[i] = GenericRadialLookup::soften(
        (p_.radial_nominal_stiffness > .5
           ? (i < 2 ? p_.c_t_f : p_.c_t_r)
           : radial_lookup_.stiffness(
               i >= 2, v, p_.inflation_absolute_Pa[i] - p_.atmosphere_Pa, f)) *
          p_.radial_stiffness_scale,
        p_.radial_softening_fraction, fx,
        i < 2 ? p_.radial_force_reference_f : p_.radial_force_reference_r);
      if (p_.radial_conditioning_enabled > .5 && p_.fy_radial_mode > .5) {
        const auto & t = i < 2 ? p_.vertical_front : p_.vertical_rear;
        const double fy = p_.radial_condition_fy_N[i] * (i % 2 ? -1. : 1.);
        const double rho = f / std::max(1000., radial_k_Npm_[i]);
        const double factor =
          p_.fy_radial_mode < 1.5
            ? manufacturer_fy_factor(
                t, p_.inflation_absolute_Pa[i] - p_.atmosphere_Pa, v, -p_.camber_out_rad[i], rho,
                fy, f)
            : sidewall_radial_factor(
                fy / (i < 2 ? p_.lateral_stiffness_f : p_.lateral_stiffness_r),
                i < 2 ? p_.sidewall_height_f : p_.sidewall_height_r, f, radial_k_Npm_[i]);
        radial_k_Npm_[i] *= factor;
      }
      radial_force_coefficient_Npm_[i] = radial_k_Npm_[i];
      if (p_.radial_conditioning_enabled < .5) {
        radial_compliance_[i].reset(
          radial_lookup_, i >= 2, v, p_.inflation_absolute_Pa[i] - p_.atmosphere_Pa);
        const double factor = GenericRadialLookup::soften(
          1., p_.radial_softening_fraction, fx,
          i < 2 ? p_.radial_force_reference_f : p_.radial_force_reference_r);
        radial_force_coefficient_Npm_[i] =
          f > 1e-8
            ? (p_.radial_nominal_stiffness > .5
                 ? (i < 2 ? p_.c_t_f : p_.c_t_r) * factor * p_.radial_stiffness_scale
                 : f * factor / radial_compliance_[i].values(f)[0] * p_.radial_stiffness_scale)
            : radial_k_Npm_[i];
      }
    }
  }
  ScopedStageRoadCache road_cache(road_surface_, p_.radial_robust_solver_enabled > .5);
  radial_iteration_residual_ = 0.;
  radial_root_fallback_used_ = 0.;
  radial_root_evaluations_ = 0.;
  auto evaluate_radial_contact = [&]() {
    calculate_tire_spring_initial_compression_m();
    calculate_tire_spring_force_N();
#ifdef EAV24_STEP_TIMING
    const auto contact_start = step_timing_enabled_ ? TimingClock::now() : TimingClock::time_point{};
#endif
    calculate_local_road_contact();
#ifdef EAV24_STEP_TIMING
    if (step_timing_enabled_) {
      step_timing_.contact_us += elapsed_us(contact_start);
      ++step_timing_.contact_calls;
    }
#endif
    calculate_vertical_tire_force_N();
    if (
      (p_.mf62_vertical_enabled > .5 && p_.mf62_rolling_enabled > .5) ||
      p_.effective_radius_only > .5)
      for (int i = 0; i < 4; ++i) {
        const auto & t = i < 2 ? p_.vertical_front : p_.vertical_rear;
        imr_.dynamic_tire_radius_m[i] = mf62_effective_radius(
          t, imr_.vertical_tire_force_N[i], p_.inflation_absolute_Pa[i] - p_.atmosphere_Pa,
          wheel_speeds_radps_[i]);
      }
    calculate_velocity_wheel_over_ground_mps();
    calculate_velocity_wheel_over_ground_tire_frame_mps();
    calculate_velocity_tire_rotation_mps();
    calculate_tire_longitudinal_slip();
    calculate_tire_slip_angle_rad();
    calculate_tire_forces_tire_frame_N();
  };
  auto initial_coefficients = radial_force_coefficient_Npm_;
  auto initial_tangents = radial_k_Npm_;
  auto legacy_iteration = [&]() {
    for (int iteration = 0; iteration < 40; ++iteration) {
      evaluate_radial_contact();
      if (p_.radial_lookup_enabled < .5 || p_.radial_conditioning_enabled > .5) break;
      double maximum = 0.;
      double_per_wheel_t next{};
      for (int i = 0; i < 4; ++i) {
        const double speed = radial_lookup_speed_mps();
        const double force = std::max(0., imr_.vertical_tire_force_N[i]);
        const double k =
          (p_.radial_nominal_stiffness > .5
             ? (i < 2 ? p_.c_t_f : p_.c_t_r)
             : radial_lookup_.stiffness(
                 i >= 2, speed, p_.inflation_absolute_Pa[i] - p_.atmosphere_Pa, force)) *
          p_.radial_stiffness_scale;
        double factor = GenericRadialLookup::soften(
          1., p_.radial_softening_fraction, imr_.tire_forces_tire_frame_N[i][0],
          i < 2 ? p_.radial_force_reference_f : p_.radial_force_reference_r);
        if (p_.fy_radial_mode > .5) {
          const double fy = imr_.tire_forces_tire_frame_N[i][1] * (i % 2 ? -1. : 1.);
          const double rho = p_.radial_nominal_stiffness > .5
                               ? force / std::max(1000., k)
                               : radial_compliance_[i].values(force)[0];
          const double change =
            p_.fy_radial_mode < 1.5
              ? manufacturer_fy_factor(
                  i < 2 ? p_.vertical_front : p_.vertical_rear,
                  p_.inflation_absolute_Pa[i] - p_.atmosphere_Pa, speed,
                  -(p_.camber_out_rad[i] + kinematic_camber_rad_[i]), rho, fy, force)
              : sidewall_radial_factor(
                  x_vec_[x::lateral_deflection_tire_FL_m + i],
                  i < 2 ? p_.sidewall_height_f : p_.sidewall_height_r, force, k);
          if (!(change > .05 && std::isfinite(change)))
            throw std::domain_error("Fy radial sensitivity leaves positive domain");
          factor *= change;
        }
        radial_k_Npm_[i] = k * factor;
        next[i] = force > 1e-8
                    ? (p_.radial_nominal_stiffness > .5
                         ? (i < 2 ? p_.c_t_f : p_.c_t_r) * factor * p_.radial_stiffness_scale
                         : force * factor / radial_compliance_[i].values(force)[0] *
                             p_.radial_stiffness_scale)
                    : radial_k_Npm_[i];
        maximum = std::max(maximum, std::abs(next[i] / radial_force_coefficient_Npm_[i] - 1.));
      }
      radial_iteration_residual_ = maximum;
      if (maximum < 1e-9) break;
      if (iteration == 39)
        throw std::runtime_error("Radial constitutive iteration did not converge");
      for (int i = 0; i < 4; ++i)
        radial_force_coefficient_Npm_[i] = .5 * (radial_force_coefficient_Npm_[i] + next[i]);
    }
  };
  try {
    legacy_iteration();
  } catch (const std::exception & error) {
    const std::string reason = error.what();
    const bool eligible =
      reason.find("Radial constitutive iteration") != std::string::npos ||
      reason.find("Radial softening outside positive-force domain") != std::string::npos ||
      reason.find("Fy radial sensitivity") != std::string::npos ||
      reason.find("Manufacturer Fy correction") != std::string::npos;
    if (!(p_.radial_robust_solver_enabled > .5 && p_.radial_lookup_enabled > .5 &&
          p_.radial_conditioning_enabled < .5 && eligible))
      throw;
    radial_root_fallback_used_ = 1.;
    auto failed_coefficients = radial_force_coefficient_Npm_;
    radial_force_coefficient_Npm_ = initial_coefficients;
    radial_k_Npm_ = initial_tangents;
    auto target = [&](int i) {
      const double force = std::max(0., imr_.vertical_tire_force_N[i]);
      const double speed = radial_lookup_speed_mps();
      const double base =
        (p_.radial_nominal_stiffness > .5
           ? (i < 2 ? p_.c_t_f : p_.c_t_r)
           : radial_lookup_.stiffness(
               i >= 2, speed, p_.inflation_absolute_Pa[i] - p_.atmosphere_Pa, force)) *
        p_.radial_stiffness_scale;
      const double rho = p_.radial_nominal_stiffness > .5 ? force / std::max(1000., base)
                                                          : radial_compliance_[i].values(force)[0];
      const double coefficient =
        p_.radial_nominal_stiffness > .5
          ? base
          : (force > 1e-8 ? force / rho * p_.radial_stiffness_scale : base);
      const double fx = imr_.tire_forces_tire_frame_N[i][0],
                   reference = i < 2 ? p_.radial_force_reference_f : p_.radial_force_reference_r;
      const double fxfactor = 1. - p_.radial_softening_fraction * std::pow(fx / reference, 2);
      double fyfactor = 1.;
      if (p_.fy_radial_mode > .5) {
        const double fy = imr_.tire_forces_tire_frame_N[i][1] * (i % 2 ? -1. : 1.);
        fyfactor = p_.fy_radial_mode < 1.5
                     ? manufacturer_fy_factor(
                         i < 2 ? p_.vertical_front : p_.vertical_rear,
                         p_.inflation_absolute_Pa[i] - p_.atmosphere_Pa, speed,
                         -(p_.camber_out_rad[i] + kinematic_camber_rad_[i]), rho, fy, force, false)
                     : sidewall_radial_factor(
                         x_vec_[x::lateral_deflection_tire_FL_m + i],
                         i < 2 ? p_.sidewall_height_f : p_.sidewall_height_r, force, base);
      }
      if (!(std::isfinite(fxfactor) && std::isfinite(fyfactor)))
        throw std::domain_error("Nonfinite radial root residual");
      return std::array<double, 4>{
        coefficient * fxfactor * fyfactor, fxfactor, fyfactor, base * fxfactor * fyfactor};
    };
    for (int i = 0; i < 4; ++i) {
      const double base = std::max(initial_coefficients[i], initial_tangents[i]);
      const double minimum = base * 1e-4;
      double maximum = base * (p_.fy_radial_mode > .5 && p_.fy_radial_mode < 1.5 ? 5. : 1.001) *
                       (p_.radial_nominal_stiffness > .5 ? 1. : 10.);
      if (p_.belt_enabled > .5) {
        const double mass = i < 2 ? p_.belt_mass_f : p_.belt_mass_r,
                     frequency = i < 2 ? p_.belt_frequency_f : p_.belt_frequency_r;
        maximum = std::min(maximum, mass * std::pow(2. * M_PI * frequency, 2) * .999);
      }
      const auto coupled_seed_coefficients=radial_force_coefficient_Npm_;
      const auto coupled_seed_tangents=radial_k_Npm_;
      auto coupled_seed_radii=imr_.dynamic_tire_radius_m;
      auto residual = [&](double coefficient) {
        radial_force_coefficient_Npm_[i] = coefficient;
        // The lookup's secant force coefficient and tangent are distinct.
        // Belt damping conversion uses the tangent, so equilibrate that small
        // algebraic feedback at each bracket trial rather than setting Ktan=Csec.
        if (p_.radial_nominal_stiffness > .5)
          radial_k_Npm_[i] = coefficient;
        else
          radial_k_Npm_[i] = initial_tangents[i];
        for (int tangent_iteration = 0; tangent_iteration < 60; ++tangent_iteration) {
          ++radial_root_evaluations_;
          evaluate_radial_contact();
          auto t = target(i);
          if (p_.radial_nominal_stiffness > .5 || p_.belt_enabled < .5 || t[3] <= 0.)
            return coefficient - t[0];
          if (std::abs(radial_k_Npm_[i] / t[3] - 1.) < 1e-12) return coefficient - t[0];
          radial_k_Npm_[i] = .5 * (radial_k_Npm_[i] + t[3]);
        }
        // A relaxed fixed point need not be contractive. Solve its original
        // tangent equation on the positive series-compliance domain instead.
        const double radius_seed=imr_.dynamic_tire_radius_m[i];
        auto tangent_residual=[&](double tangent){
          radial_k_Npm_[i]=tangent;
          imr_.dynamic_tire_radius_m[i]=radius_seed;
          // Ring camber shifts contact queries by radius*gamma. Close this
          // weak effective-radius feedback too, so each scalar trial has a
          // repeatable constitutive value rather than a previous-trial radius.
          for(int geometry_iteration=0;geometry_iteration<24;++geometry_iteration){
            const double before=imr_.dynamic_tire_radius_m[i];
            ++radial_root_evaluations_;evaluate_radial_contact();
            if(std::abs(imr_.dynamic_tire_radius_m[i]-before)<1e-12*std::max(.1,std::abs(before)))
              return tangent-target(i)[3];
          }
          // A rolling-radius/contact fixed point can also lose contraction.
          // Solve its unchanged geometric equation on the positive radius
          // domain, bounded above by the unloaded MF radius at this spin.
          const auto &vertical=i<2?p_.vertical_front:p_.vertical_rear;
          const double free_radius=mf62_free_radius(vertical,wheel_speeds_radps_[i]);
          auto radius_residual=[&](double radius){
            imr_.dynamic_tire_radius_m[i]=radius;
            ++radial_root_evaluations_;evaluate_radial_contact();
            return radius-imr_.dynamic_tire_radius_m[i];
          };
          try{
            const double radius=positive_residual_root(
              radius_residual,radius_seed,1e-6,free_radius*(1.+1e-12));
            const double error=radius_residual(radius);
            if(std::abs(error)>1e-12*std::max(.1,std::abs(radius)))
              throw std::runtime_error("Rolling-radius root fails geometric residual");
            return tangent-target(i)[3];
          }catch(const std::exception &failure){
            throw std::runtime_error(std::string("Radial TOMS748: radius feedback root failed: ")+failure.what());
          }
        };
        try{
          const double tangent=positive_residual_root(tangent_residual,initial_tangents[i],minimum,maximum);
          tangent_residual(tangent);
          return coefficient-target(i)[0];
        }catch(const std::exception &failure){
          throw std::runtime_error(std::string("Radial TOMS748: tangent feedback root failed: ")+failure.what());
        }
      };
      try{
        const double seed = std::clamp(failed_coefficients[i], minimum, maximum);
        double lo = std::max(minimum, seed * .8), hi = std::min(maximum, seed * 1.2);
        double fl = residual(lo), fh = residual(hi);
        for (int expansion = 0; fl * fh > 0. && expansion < 20; ++expansion) {
          const double a = std::max(minimum, lo * .5), b = std::min(maximum, hi * 2.);
          if (a == lo && b == hi) break;
          lo = a;
          hi = b;
          fl = residual(lo);
          fh = residual(hi);
        }
        if (!(lo < hi && std::isfinite(fl) && std::isfinite(fh) && fl * fh <= 0.))
          throw std::domain_error("Radial TOMS748: no sign-changing positive coefficient bracket");
        boost::uintmax_t iterations = 80;
        auto bracket = boost::math::tools::toms748_solve(
          residual, lo, hi, fl, fh, boost::math::tools::eps_tolerance<double>(44), iterations);
        const double root = .5 * (bracket.first + bracket.second);
        const double rr = residual(root);
        auto t = target(i);
        if (!(std::abs(rr) / std::max(1., root) < 1e-9 && t[1] > .05 && t[2] > .05 && t[2] < 5.))
          throw std::domain_error(
            "Radial TOMS748: root fails residual or unchanged positive-force guards");
        radial_force_coefficient_Npm_[i] = root;
        radial_k_Npm_[i] = t[3];
      }catch(const std::exception &scalar_failure){
        const std::string reason=scalar_failure.what();
        if(reason.find("Radial TOMS748")==std::string::npos&&
           reason.find("Positive scalar")==std::string::npos&&
           reason.find("Tire collider leaves measured surface coverage")==std::string::npos)throw;
        // Nested scalar trials may have no admissible intermediate radius even
        // when the original three-equation system has a local positive root.
        // Reset every trial completely and solve the unchanged equations jointly.
        const double cscale=std::max(1000.,initial_coefficients[i]);
        const double kscale=std::max(1000.,initial_tangents[i]);
        const auto &vertical=i<2?p_.vertical_front:p_.vertical_rear;
        const double rscale=mf62_free_radius(vertical,wheel_speeds_radps_[i]);
        Eigen::Vector3d lower(minimum/cscale,minimum/kscale,1e-6/rscale);
        Eigen::Vector3d upper(maximum/cscale,maximum/kscale,1.+1e-12);
        Eigen::Vector3d seed(initial_coefficients[i]/cscale,initial_tangents[i]/kscale,
                             coupled_seed_radii[i]/rscale);
        auto coupled=[&](const Eigen::Vector3d &u){
          radial_force_coefficient_Npm_=coupled_seed_coefficients;
          radial_k_Npm_=coupled_seed_tangents;
          imr_.dynamic_tire_radius_m=coupled_seed_radii;
          const double coefficient=u[0]*cscale,tangent=u[1]*kscale,radius=u[2]*rscale;
          radial_force_coefficient_Npm_[i]=coefficient;
          radial_k_Npm_[i]=tangent;
          imr_.dynamic_tire_radius_m[i]=radius;
          ++radial_root_evaluations_;evaluate_radial_contact();
          const auto desired=target(i);
          return Eigen::Vector3d((coefficient-desired[0])/std::max(1.,coefficient),
            (tangent-desired[3])/std::max(1.,tangent),
            (radius-imr_.dynamic_tire_radius_m[i])/(1e-3*std::max(.1,std::abs(radius))));
        };
        try{
          Eigen::Vector3d root;
          try{root=coupled_radial_root(coupled,seed,lower,upper);}
          catch(const std::exception&){root=hybrid_radial_root(coupled,seed,lower,upper);}
          const Eigen::Vector3d error=coupled(root);const auto desired=target(i);
          if(!(error.lpNorm<Eigen::Infinity>()<1e-9))
            throw std::runtime_error("Final coupled residual rejected");
          if(!(desired[1]>.05&&desired[2]>.05&&desired[2]<5.))
            throw std::domain_error("Radial softening outside positive-force domain after coupled solve");
        }catch(const std::domain_error &failure){
          if(std::string(failure.what()).find("Radial softening outside positive-force domain")!=std::string::npos)throw;
          throw std::runtime_error(std::string("Radial TOMS748: coupled root failed: ")+failure.what());
        }catch(const std::exception &failure){
          throw std::runtime_error(std::string("Radial TOMS748: coupled root failed: ")+failure.what());
        }
      }
    }
    evaluate_radial_contact();
    radial_iteration_residual_ = 0.;
    for (int i = 0; i < 4; ++i) {
      auto t = target(i);
      radial_iteration_residual_ = std::max(
        radial_iteration_residual_, std::abs(t[0] / radial_force_coefficient_Npm_[i] - 1.));
      if (!(t[1] > .05 && t[2] > .05 && t[2] < 5.))
        throw std::domain_error("Radial TOMS748: final force law outside unchanged guards");
    }
    if (radial_iteration_residual_ >= 1e-9)
      throw std::runtime_error("Radial TOMS748: coupled final residual exceeds tolerance");
  }

  calculate_tire_forces_N();
  if(p_.ring_enabled>.5)calculate_tire_rolling_resistance_N();
  calculate_ring_reactions();
  for (int i = 0; i < 4; ++i) {
    if (p_.ring_enabled < .5) {
      imr_.tire_forces_N[i][0] += imr_.contact_force_x_N[i];
      imr_.tire_forces_N[i][1] += imr_.contact_force_y_N[i];
    } // The ring branch has already transmitted these loads through its sidewall.
    const double delta = imr_.effective_steering_angle_per_wheel_rad[i];
    const double lateral =
      -std::sin(delta) * imr_.tire_forces_N[i][0] + std::cos(delta) * imr_.tire_forces_N[i][1];
    imr_.tire_lateral_target_m[i] =
      lateral / (i < 2 ? p_.lateral_stiffness_f : p_.lateral_stiffness_r);
  }
  calculate_suspension_spring_initial_compression_m();
  calculate_suspension_spring_compression_m();
  calculate_suspension_damper_compression_speed_mps();
  calculate_suspension_spring_force_N();
  for (int i = 0; i < 4; ++i) {
    const double compression = std::max(
      0.,
      -(imr_.suspension_spring_compression_m[i] + imr_.suspension_spring_initial_compression_m[i]) -
        (i < 2 ? p_.bump_gap_f : p_.bump_gap_r));
    const double actual =
      p_.kinematics_enabled > .5
        ? std::max(
            0., -geometry_equivalent_extension_m_[i] - (i < 2 ? p_.bump_gap_f : p_.bump_gap_r))
        : compression;
    imr_.suspension_bump_force_N[i] = -(i < 2 ? p_.bump_stiffness_f : p_.bump_stiffness_r) *
                                      actual * actual *
                                      (p_.kinematics_enabled > .5 ? geometry_leverage_[i] : 1.);
    imr_.endstop_force_N[i] = 0.;
    if (p_.endstop_enabled > .5) {
      const double ratio = i < 2 ? p_.damper_motion_ratio_f : p_.damper_motion_ratio_r;
      const double bump = i < 2 ? p_.stroke_bump_front : p_.stroke_bump_rear;
      const double rebound = i < 2 ? p_.stroke_rebound_front : p_.stroke_rebound_rear;
      if (!(ratio > 0. && bump > 0. && rebound > 0. && p_.endstop_stiffness > 0. &&
            p_.endstop_damping >= 0.))
        throw std::invalid_argument("Invalid suspension travel limits");
      const double q = geometry_equivalent_extension_m_[i] / ratio;
      const double v =
        imr_.suspension_damper_compression_speed_mps[i] * geometry_leverage_[i] / ratio;
      const double dx = std::max(0., q - rebound) - std::max(0., -q - bump);
      double force = p_.endstop_stiffness * dx;
      if ((q > rebound && v > 0.) || (q < -bump && v < 0.)) force += p_.endstop_damping * v;
      imr_.endstop_force_N[i] = force * geometry_leverage_[i] / ratio;
      imr_.suspension_bump_force_N[i] += imr_.endstop_force_N[i];
    }
  }
  // One third element per axle: its compression is the MEAN wheel jounce
  // divided by wheel/device ratio. Each wheel receives half the device force
  // divided by that same ratio (virtual work). Subtract the nominal preload
  // already carried by the native static spring equilibrium.
  for (int axle = 0; axle < 2; ++axle) {
    const int left = 2 * axle, right = left + 1;
    const double stiffness = axle == 0 ? p_.heave_stiffness_f : p_.heave_stiffness_r;
    const double coil = axle == 0 ? p_.heave_spring_f : p_.heave_spring_r;
    const double motion_ratio = axle == 0 ? p_.heave_motion_ratio_f : p_.heave_motion_ratio_r;
    const double gap = axle == 0 ? p_.heave_gap_f : p_.heave_gap_r;
    if (!(std::isfinite(stiffness) && stiffness >= 0. && std::isfinite(coil) && coil >= 0. &&
          std::isfinite(motion_ratio) && motion_ratio > 0. && std::isfinite(gap)))
      throw std::invalid_argument("Third-element stiffness/ratio/gap invalid");
    if (stiffness == 0. && coil == 0.) {
      imr_.suspension_heave_force_N[left] = imr_.suspension_heave_force_N[right] = 0.;
      imr_.suspension_heave_energy_J[left] = imr_.suspension_heave_energy_J[right] = 0.;
      continue;
    }
    const double jounce =
      p_.kinematics_enabled > .5
        ? (axle == 0 ? .5 * (kinematic_spring_m_[left] + kinematic_spring_m_[right])
                     : -.5 * (geometry_equivalent_extension_m_[left] +
                              geometry_equivalent_extension_m_[right]))
        : -.5 * (imr_.suspension_spring_compression_m[left] +
                 imr_.suspension_spring_initial_compression_m[left] +
                 imr_.suspension_spring_compression_m[right] +
                 imr_.suspension_spring_initial_compression_m[right]);
    double device_travel = jounce / motion_ratio;
    double jl = .5 / motion_ratio, jr = jl;
    if (p_.kinematics_enabled > .5) {
      if (axle == 0) {
        jl *= kinematic_spring_jac_[left];
        jr *= kinematic_spring_jac_[right];
      } else {
        const auto slider = suspension_smooth_lookup(
          .5 * (kinematic_arb_m_[left] + kinematic_arb_m_[right]), p_.rear_third_bar_travel,
          p_.rear_third_device_compression);
        device_travel = slider.first;
        jl = .5 * slider.second * kinematic_arb_jac_[left];
        jr = .5 * slider.second * kinematic_arb_jac_[right];
      }
    }
    const double compression = std::max(0., device_travel - gap), initial = std::max(0., -gap);
    double force =
      -(stiffness * (compression * compression - initial * initial) + coil * device_travel);
    double energy = stiffness * ((std::pow(compression, 3) - std::pow(initial, 3)) / 3. -
                                 initial * initial * device_travel) +
                    .5 * coil * device_travel * device_travel;
    if (p_.endstop_enabled > .5 && p_.kinematics_enabled > .5) {
      const double initial_stroke = axle == 0 ? p_.endstop_heave_initial_stroke_f : p_.endstop_heave_initial_stroke_r,
                   travel = axle == 0 ? p_.endstop_heave_travel_f : p_.endstop_heave_travel_r,
                   c = initial_stroke + device_travel;
      const double penetration = std::max(0., c - travel) - std::max(0., -c);
      force -= p_.endstop_stiffness * penetration;
      energy += .5 * p_.endstop_stiffness * penetration * penetration;
    }
    imr_.suspension_heave_force_N[left] = force * jl;
    imr_.suspension_heave_force_N[right] = force * jr;
    imr_.suspension_heave_energy_J[left] = imr_.suspension_heave_energy_J[right] = .5 * energy;
  }
  calculate_suspension_damper_force_N();
  calculate_antiroll_bar_force_N();
  calculate_tire_rolling_resistance_N();
  calculate_axle_vertical_force_N();
  calculate_resulting_suspension_force_N();
  calculate_resulting_vertical_force_on_wheel_N();
  calculate_aerodynamics();
  calculate_resulting_force_N();
  calculate_resulting_torque_Nm();
  if(p_.ring_enabled>.5){imr_.resulting_torque_Nm.x+=ring_pitch_correction_[0];imr_.resulting_torque_Nm.y+=ring_pitch_correction_[1];}
  // The new third element is derived from sprung-body attachment kinematics.
  // Native pitch moments use the whole-vehicle CoG and a gravity offset; apply
  // only this element's moment correction so its force obeys virtual work.
  if (p_.kinematics_enabled < .5)
    imr_.resulting_torque_Nm.y +=
      (imr_.suspension_heave_force_N.front_left + imr_.suspension_heave_force_N.front_right) *
        (p_dep_.l_f_sprung - p_.l_f) -
      (imr_.suspension_heave_force_N.rear_left + imr_.suspension_heave_force_N.rear_right) *
        (p_dep_.l_r_sprung - p_dep_.l_r);
  for (int i = 0; i < 4; ++i) imr_.resulting_torque_Nm.z += imr_.contact_yaw_moment_Nm[i];
  // clang-format on
}
}  // namespace tam::ocd::vehicle_dynamics
