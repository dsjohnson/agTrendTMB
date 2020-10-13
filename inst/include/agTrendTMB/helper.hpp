
using namespace density;

// Utility functions
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
bool isFinite(Type x){
  return R_finite(asDouble(x));
}


template<class Type>
Type iou_corr(Type r, Type s, Type phi){
  Type u = sqrt((r-s)*(r-s));
  Type out = (2*phi*std::min(s,r) + exp(-phi*s) + exp(-phi*r) - exp(-phi*u) - 1)/pow(phi,2);
  return(out);
}


template<class Type>
Type gp_ll_matern(vector<Type> x, vector<Type> t, Type ln_sigma, Type ln_phi, Type order){
  int n=x.size();
  Type kappa = order + 0.5;
  Type phi = exp(ln_phi);
  matrix<Type> C(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      Type u = abs(t(i)-t(j));
      C(i,j) = exp(2*ln_sigma) * matern(u, phi, kappa);
    }
  }
  Type ll = -1.0*density::MVNORM_t<Type>(C)(x);
  return(ll);
}

template<class Type>
Type gp_ll_iou(vector<Type> x, vector<Type> t, Type ln_sigma, Type ln_phi){
  int n=x.size();
  Type phi = exp(ln_phi);
  matrix<Type> C(n,n);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      Type r = t(i);
      Type s = t(j);
      C(i,j) = exp(2*ln_sigma) * iou_corr(r, s, phi);
    }
    C(i,i) += Type(1.0E-8);
  }
  Type ll = -1.0*density::MVNORM_t<Type>(C)(x);
  return(ll);
}

template<class Type>
Type gp_ll_ibm(vector<Type> x, vector<Type> t, Type ln_sigma){
  int n=x.size();
  matrix<Type> C(n,n);
  for(int i=0;i<n;i++){
    for(int j=i;j<n;j++){
      Type r = t(i);
      Type s = t(j);
      C(i,j) = exp(2*ln_sigma) * (r*r*(s/2 - r/6));
      C(j,i) = C(i,j);
    }
    C(i,i) += Type(1.0E-8);
  }
  Type ll = -1.0*density::MVNORM_t<Type>(C)(x);
  return(ll);
}

// GRMF log-likelihood
// template<class Type>
// Type gmrf_ll(vector<Type> alpha, SparseMatrix<Type> Q, Type Qrank, Type log_lambda){
//   Type lambda = exp(log_lambda);
//   Type out = Type(0.5)*Qrank*log_lambda - 0.5*lambda*density::GMRF(Q).Quadform(alpha);
// }
// mvnorm for known Q
template<class Type>
Type mnorm_Q_ll(vector<Type> x, vector<Type> m, matrix<Type> Q){
  vector<Type> r = x-m;
  Type out = - 0.5*(r*(Q*r)).sum();
  return(out);
}


template <class Type>
Type disc_norm(Type x, Type mu, Type sd){
  Type f = 0.0;
  f = pnorm(Type(x+0.5), mu, sd)-pnorm(Type(x-0.5), mu, sd);
  if(x==0){
    f += pnorm(Type(-0.5), mu, sd);
  }
  f = log(f);
  return(f);
}

template <class Type>
Type disc_tr_norm(Type x, Type mu, Type sd, Type a, Type b){
  Type f = 0.0;
  Type K = pnorm(b, mu, sd)-pnorm(a, mu, sd);
  f = (pnorm(x+0.5, mu, sd) - pnorm(x-0.5, mu, sd))/K;
  if(x==0){
    f += (pnorm(Type(-0.5), mu, sd)-pnorm(a, mu, sd))/K;
  }
  f = log(f);
  return(f);
}

template <class Type>
Type disc_tr_zinorm(Type x, Type mu, Type sd, Type p, Type b){
  Type f = 0.0;
  Type K = pnorm(b, mu, sd)-pnorm(0.5, mu, sd);
  if(x>0){
    f = (1-p)*(pnorm(x+0.5, mu, sd) - pnorm(x-0.5, mu, sd))/K;
  } else if(x==0){
    f = p;
  }
  f = log(f);
  return(f);
}

template <class Type>
Type disc_mix_norm(Type x, vector<Type> mu, vector<Type> sd, vector<Type> w, vector<Type> delta, vector<int> use_zi){
  int K = mu.size();
  vector<Type> f(K); f.SetZero();
  vector<Type> zi(K); zi.SetZero();
  for(int i=0; i<K; i++){
    if(use_zi(i)){
      zi(i) = pnorm(Type(-0.5-delta(i)), mu(i), sd(i));
    }
    f(i) = (pnorm(Type(x+0.5), mu(i), sd(i)) - pnorm(Type(x-0.5), mu(i), sd(i)))/(1-zi(i));
    if(x==0){
      f(i) += (pnorm(Type(-0.5), mu(i), sd(i)) - zi(i))/(1-zi(i));
    }
    f = log(w(i)) + log(f);
  }
  return(f.sum());
}


template <class Type>
Type disc_lognorm(Type x, Type mu, Type sd, Type a, Type b){
  Type f = 0.0;
  if(a>log(0.5)) return(-Type(INFINITY));
  if(b<log(x+0.5)) return(-Type(INFINITY));
  Type K = pnorm(b, mu, sd)-pnorm(a,mu,sd);
  if(x>0){
    f = (pnorm(Type(log(x+0.5)), mu, sd) - pnorm(Type(log(x-0.5)), mu, sd))/K;
  }
  else if(x==0){
    f = (pnorm(Type(log(0.5)), mu, sd) - pnorm(a, mu, sd))/K;
  } 
  f = log(f);
  return(f);
}

template <class Type>
Type disc_zilognorm(Type x, Type mu, Type sd, Type p, Type b){
  Type f = 0.0;
  if(b<log(x+0.5)) return(Type(INFINITY));
  Type K = pnorm(b, mu, sd)-pnorm(Type(0.5), mu, sd);
  if(x>0){
    f = (1-p)*(pnorm(Type(log(x+0.5)), mu, sd) - pnorm(Type(log(x-0.5)), mu, sd))/K;
  }
  else if(x==0){f = p;} 
  f = log(f);
  return(f);
}


template<class Type>
Type log_growth_prior(Type x, Type ul){
  Type out;
  if(x>=0){
    out = dnorm(x, Type(0), Type(ul/2), true);
  } else {
    out = dnorm(Type(0), Type(0), Type(ul/2), true);
  }
  return(out);
}
