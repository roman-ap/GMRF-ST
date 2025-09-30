data {
  int<lower=1> N;
  int<lower=1> P;
  array[N] int<lower=0> y;
  vector[N] log_offset;
  matrix[N, P] X;
  matrix[N, N] A;
  matrix[N, N] D;
  vector[N] lambda;
}

transformed data {
  vector[N] zeros;
  zeros = rep_vector(0, N);
}

parameters {
  real alpha;
  vector[P] beta;
  vector[N] phi;
  real<lower=0> phi_tau;
  real<lower=1/min(lambda), upper=1/max(lambda)> phi_rho;
}

transformed parameters {
  real<lower=0> phi_sigma;
  phi_sigma = inv(sqrt(phi_tau));
}

model {
  target += poisson_log_lpmf(y | log_offset + phi);
  target += multi_normal_prec_lpdf(phi | alpha + X * beta, phi_tau * (D - phi_rho * A));
  target += normal_lpdf(alpha | 1, 1);
  target += normal_lpdf(beta | 2, 1);
  target += gamma_lpdf(phi_tau | 1000, 10);
  target += uniform_lpdf( phi_rho | 0, 1);
}
