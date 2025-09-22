functions {
#include ST_lpdf.stan
}

data {
  int<lower=1> S;
  array[S] int<lower=0> y;
  vector[S] log_offset;
  int<lower=1> S_sparse;
  vector[S_sparse] A_w;
  array[S_sparse] int A_v;
  array[S+1] int A_u;
  vector[S] D_ii;
  real log_det_D;
  vector[S] lambda;
}

transformed data {
  vector[S] zeros = rep_vector(0, S);
}

parameters {
  vector[S] phi;
  real<lower=1/min(lambda), upper=1/max(lambda)> phi_rho;
}

model {
  target += poisson_log_lpmf(y | log_offset + phi);
  target += wcar_normal_lpdf(phi | zeros, phi_rho, A_w, A_v, A_u, D_ii, log_det_D, lambda, S);
}
