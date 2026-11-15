// Copyright 2026 TUM Autonomous Motorsport
#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <boost/math/tools/toms748_solve.hpp>
namespace tam::ocd::vehicle_dynamics {
// Search a positive constitutive coordinate; do not clip its residual or accept
// a narrow bracket across a jump without checking the original equation.
template<class Function>
double positive_residual_root(Function residual,double seed,double minimum,double maximum) {
  if(!(std::isfinite(seed)&&minimum>0.&&minimum<maximum&&std::isfinite(maximum)))
    throw std::invalid_argument("Invalid positive scalar root domain");
  seed=std::clamp(seed,minimum,maximum);
  double lo=std::max(minimum,.8*seed),hi=std::min(maximum,1.2*seed);
  double fl=residual(lo),fh=residual(hi);
  for(int expansion=0;fl*fh>0.&&expansion<24;++expansion){
    const double a=std::max(minimum,.5*lo),b=std::min(maximum,2.*hi);
    if(a==lo&&b==hi)break;
    lo=a;hi=b;fl=residual(lo);fh=residual(hi);
  }
  if(!(lo<hi&&std::isfinite(fl)&&std::isfinite(fh)&&fl*fh<=0.))
    throw std::runtime_error("No sign-changing positive scalar bracket");
  boost::uintmax_t iterations=100;
  const auto bracket=boost::math::tools::toms748_solve(residual,lo,hi,fl,fh,boost::math::tools::eps_tolerance<double>(44),iterations);
  const double root=.5*(bracket.first+bracket.second),error=residual(root);
  if(!(std::isfinite(error)&&std::abs(error)/std::max(1.,std::abs(root))<1e-9))
    throw std::runtime_error("Positive scalar bracket closed without a valid residual");
  return root;
}
}
