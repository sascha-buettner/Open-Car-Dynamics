// Copyright 2026 TUM Autonomous Motorsport
#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/QR>
namespace tam::ocd::vehicle_dynamics {
// Dimensionless, bounded simultaneous solve. Bounds restrict trial coordinates;
// no force law or residual is clipped. Failure never returns an approximate root.
template<class Function>
Eigen::Vector3d coupled_radial_root(Function residual,Eigen::Vector3d seed,
  const Eigen::Vector3d &lower,const Eigen::Vector3d &upper) {
  if(!(lower.allFinite()&&upper.allFinite()&&(lower.array()>0.).all()&&
       (lower.array()<upper.array()).all()&&seed.allFinite()))
    throw std::invalid_argument("Invalid coupled radial root domain");
  auto bound=[&](Eigen::Vector3d v){return v.cwiseMax(lower).cwiseMin(upper).eval();};
  auto evaluate=[&](const Eigen::Vector3d &v,Eigen::Vector3d &f){
    try{f=residual(v);return f.allFinite();}catch(const std::exception&){return false;}
  };
  Eigen::Vector3d x=bound(seed),f;
  if(!evaluate(x,f))throw std::runtime_error("Coupled radial seed is outside evaluable domain");
  for(int iteration=0;iteration<60;++iteration){
    if(f.lpNorm<Eigen::Infinity>()<3e-10)return x;
    Eigen::Matrix3d jacobian;
    for(int column=0;column<3;++column){
      const double h=1e-6*std::max(1.,std::abs(x[column]));
      auto a=x,b=x;a[column]=std::max(lower[column],x[column]-h);
      b[column]=std::min(upper[column],x[column]+h);
      Eigen::Vector3d fa,fb;
      const bool va=a[column]<x[column]&&evaluate(a,fa);
      const bool vb=b[column]>x[column]&&evaluate(b,fb);
      if(va&&vb)jacobian.col(column)=(fb-fa)/(b[column]-a[column]);
      else if(vb)jacobian.col(column)=(fb-f)/(b[column]-x[column]);
      else if(va)jacobian.col(column)=(f-fa)/(x[column]-a[column]);
      else throw std::runtime_error("Coupled radial Jacobian has no admissible difference");
    }
    const Eigen::Matrix3d normal=jacobian.transpose()*jacobian;
    const Eigen::Vector3d gradient=jacobian.transpose()*f;
    Eigen::Vector3d diagonal=normal.diagonal().cwiseMax(1e-8);
    bool accepted=false;
    for(double ridge:{0.,1e-8,1e-6,1e-4,1e-2,1.}){
      Eigen::Vector3d step;
      if(ridge==0.)step=jacobian.colPivHouseholderQr().solve(-f);
      else{Eigen::Matrix3d matrix=normal;matrix.diagonal()+=ridge*diagonal;step=matrix.ldlt().solve(-gradient);}
      if(!step.allFinite())continue;
      const double length=step.lpNorm<Eigen::Infinity>();
      if(length>.5)step*=.5/length;
      for(double alpha=1.;alpha>=1./2048.;alpha*=.5){
        const Eigen::Vector3d candidate=bound(x+alpha*step);Eigen::Vector3d next;
        if(!evaluate(candidate,next))continue;
        if(next.squaredNorm()<f.squaredNorm()*(1.-1e-4*alpha)){
          x=candidate;f=next;accepted=true;break;
        }
      }
      if(accepted)break;
    }
    if(!accepted)throw std::runtime_error("Coupled radial residual cannot be reduced");
  }
  throw std::runtime_error("Coupled radial residual did not converge");
}
}
