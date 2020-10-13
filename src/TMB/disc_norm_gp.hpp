/// @file disc_norm_gp.hpp

#ifndef disc_norm_gp_hpp
#define disc_norm_gp_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include "agTrendTMB//helper.hpp"


/// Negative log-likelihood of the normal distribution.
template<class Type>
Type disc_norm_gp(objective_function<Type>* obj) {
  
  //Data and fixed quantities
  DATA_VECTOR(y); // data vector
  DATA_MATRIX(X); // measurement error design matrix
  DATA_MATRIX(K);
  DATA_SCALAR(lambda_tau);
  DATA_SCALAR(lambda_xi);
  DATA_INTEGER(penalty);
  
  // Parameters
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(alpha);
  PARAMETER(ln_tau);
  PARAMETER_VECTOR(xi);
  
  // Derived quantities
  int Tmax = y.size();
  vector<Type> Xbeta = (X*beta).vec();
  vector<Type> Kalpha = (K*alpha).vec();
  vector<Type> mu = Xbeta + Kalpha;
  Type tau = exp(ln_tau);
  
  vector<Type> f(7); f.setZero();// = 0.0;
  vector<Type> dll(Tmax); dll.setZero();
  vector<Type> sigma(Tmax);
  
  for(int t=0; t<Tmax; t++){
    sigma(t) = exp(xi(0) + xi(1)*mu(t));
    if(!isNA(y(t))){
      dll(t) = -disc_norm(y(t), mu(t), sigma(t));
      f(0) -= disc_norm(y(t), mu(t), sigma(t));
    }
  }
  for(int i=0; i<alpha.size(); i++){
    if(penalty){
      f(1) -= dnorm(alpha(i), Type(0), tau/sqrt(i+1), true);  
    } else{
      f(1) -= dnorm(alpha(i), Type(0), tau, true);
    }
    
  }
  f(2) -= dexp(exp(ln_tau), lambda_tau, true) + ln_tau;
  f(3) -= -xi(0);
  f(4) -= dnorm(xi(1), Type(0), lambda_xi);
  
  //Export results
  ADREPORT(mu);
  ADREPORT(xi);
  ADREPORT(ln_tau);
  REPORT(f);
  REPORT(mu);
  REPORT(alpha);
  REPORT(beta);
  REPORT(sigma);
  REPORT(xi);
  REPORT(tau);
  REPORT(dll);
  
  Type out=f.sum();
  return(out);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
