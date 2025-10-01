functions {
#include ST_lpdf.stan
}

data {
  int<lower=1> S;
  array[S] int<lower=0> y;
  vector[S] lambdas;
  real pm_mu;
  real<lower=0> psd_mu;
  real<lower=0> pm_sigma;
  real<lower=0> psd_sigma;
  real<lower=1/min(lambdas), upper=1/max(lambdas)> pm_rho;
  real<lower=0> psd_rho;
  vector[S] Ds_ii;
  real log_det_D;
  int<lower=1> As_sparse;
  vector[As_sparse] A_w;
  array[As_sparse] int A_v;
  array[S+1] int A_u;
}

transformed data {
  vector[S] zeros = rep_vector(0, S);
}

parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=1/min(lambdas), upper=1/max(lambdas)> rho;
  vector[S] phi;
}

transformed parameters {
  vector[S] eta;
  eta = mu + sigma*phi;
}

model {
  target += poisson_log_lpmf(y | eta);
  target += std_wcar_lpdf(phi | zeros, 
                                rho, 
                                log_det_D, lambdas,
                                Ds_ii,
                                A_w, A_v, A_u,
                                S);
  target += normal_lpdf(rho | pm_rho, psd_rho);
  target += normal_lpdf(mu | pm_mu, psd_mu);
  target += normal_lpdf(sigma | pm_sigma, psd_sigma);
}
