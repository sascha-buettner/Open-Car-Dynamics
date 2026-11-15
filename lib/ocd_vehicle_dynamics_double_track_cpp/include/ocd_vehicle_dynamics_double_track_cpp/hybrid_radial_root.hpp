// Copyright 2026 TUM Autonomous Motorsport
#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>
namespace tam::ocd::vehicle_dynamics {
// MINPACK-style Powell hybrid on bounded transformed coordinates. Solver status
// or a small step alone never certifies a root: check the original residual.
template<class Function>
Eigen::Vector3d hybrid_radial_root(Function residual,Eigen::Vector3d seed,
  const Eigen::Vector3d &lower,const Eigen::Vector3d &upper){
 if(!(seed.allFinite()&&lower.allFinite()&&upper.allFinite()&&
      (lower.array()>0.).all()&&(lower.array()<upper.array()).all()))
   throw std::invalid_argument("Invalid hybrid radial domain");
 struct Functor{
  Function &callback;Eigen::Vector3d lo,hi;bool balance;
  Eigen::Vector3d physical(const Eigen::VectorXd &z)const{
   Eigen::Vector3d x;
   for(int i=0;i<3;++i){const double e=std::exp(-std::abs(z[i]));
    const double p=z[i]>=0.?1./(1.+e):e/(1.+e);x[i]=lo[i]+(hi[i]-lo[i])*p;}
   return x;
  }
  int operator()(const Eigen::VectorXd &z,Eigen::VectorXd &f)const{
   try{f=callback(physical(z));if(balance)f[2]*=1e-3;return f.allFinite()?0:-1;}
   catch(const std::exception&){return -1;}
  }
 };
 Eigen::Vector3d current=seed;
 for(bool balance:{false,true}){
  Functor functor{residual,lower,upper,balance};Eigen::VectorXd z(3);
  for(int i=0;i<3;++i){const double p=std::clamp((current[i]-lower[i])/(upper[i]-lower[i]),1e-12,1.-1e-12);z[i]=std::log(p/(1.-p));}
  Eigen::HybridNonLinearSolver<Functor> solver(functor);
  solver.parameters.maxfev=600;solver.parameters.xtol=1e-12;
  solver.parameters.factor=.25;solver.parameters.epsfcn=1e-14;
  solver.solveNumericalDiff(z);current=functor.physical(z);
  try{const Eigen::Vector3d original=residual(current);
   if(original.allFinite()&&original.lpNorm<Eigen::Infinity>()<3e-10)return current;
  }catch(const std::exception&){}
 }
 throw std::runtime_error("Powell hybrid did not satisfy original coupled residual");
}
}
