functions {
#include MARCAR-functions.stan
}

data {
  int<lower=1> T;
  int<lower=1> S;
  //int<lower=1> P;
  array[T*S] int<lower=0> y;
  vector[T*S] log_offset;
  //matrix[T*S, P] X;
  real log_det_Ds;
  vector[S] lambdas;
  vector[T*S] ItDs_ii;
  int<lower=1> N_ItAs;
  vector[N_ItAs] ItAs_w;
  array[N_ItAs] int ItAs_v;
  array[T*S+1] int ItAs_u;
  int<lower=1> N_AtDs;
  vector[N_AtDs] AtDs_w;
  array[N_AtDs] int AtDs_v;
  array[T*S+1] int AtDs_u;
  int<lower=1> N_AtAs;
  vector[N_AtAs] AtAs_w;
  array[N_AtAs] int AtAs_v;
  array[T*S+1] int AtAs_u;
  vector[T*S] IItDs_ii;
  int<lower=1> N_IItAs;
  vector[N_IItAs] IItAs_w;
  array[N_IItAs] int IItAs_v;
  array[T*S+1] int IItAs_u;
}

transformed data {
  vector[T*S] zeros;
  zeros = rep_vector(0, T*S);
}

parameters {
  //vector[P] beta;
  vector[T*S] phi;
  real<lower=0> sigma_space;
  real<lower=-1, upper=1> rho_time;
  real<lower=1/min(lambdas), upper=1/max(lambdas)> rho_space;
}

model {
  target += poisson_log_lpmf(y | log_offset + phi);
  target += marcar_normal_lpdf(phi | zeros,
                               rho_time, rho_space, sigma_space,
                               log_det_Ds, lambdas,
                               ItDs_ii,
                               ItAs_w, ItAs_v, ItAs_u,
                               AtDs_w, AtDs_v, AtDs_u,
                               AtAs_w, AtAs_v, AtAs_u,
                               IItDs_ii, 
                               IItAs_w, IItAs_v, IItAs_u,
                               T, S);
  //target += normal_lpdf(beta | 1, 0.01);
  //target += normal_lpdf( rho_time | 0.9, 0.0001);
  //target += normal_lpdf( rho_space | 0.9, 0.0001);
  //target += normal_lpdf(sigma_space | 0.1, 0.01);
}

generated quantities {
  //vector[T*S] eta = X*beta + phi;
  vector[T*S] rates = exp(log_offset + phi);
}
