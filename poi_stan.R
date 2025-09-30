#####  ##### 

### Set-up
library(sf)
library(dplyr)
library(Matrix)
library(cmdstanr)
library(posterior)
library(bayesplot)

### Sweden's regions
precints <- read_sf("/proj/wildman/FixedKretsar/FixedKretsShape2023_2024.shp")
counties <- precints |> group_by(LANSKOD) |> summarise()
counties <- counties[-8,]
plot(counties)

### Simulated data for Swedish counties
# Number of counties
N<-nrow(counties)
# Binary spatial relationships  
W_match<-st_relate(counties,counties,pattern="****T1212",sparse=TRUE)
# Sparse matrices for the CAR distribution
D<-sparseMatrix(i=c( 1:N ),
                j=c( 1:N ),
                x=c( sapply(W_match,length) ),
                dims=c(N,N))
A<-sparseMatrix(i=c( rep(1:N,sapply(W_match,length)) ),
                j=c( unlist(W_match) ),
                x=c( rep(1,length(unlist(W_match))) ),
                dims=c(N,N))
C<-sparseMatrix(i=c( rep(1:N,sapply(W_match,length)) ),
                j=c( unlist(W_match) ),
                x=c( rep(1/sapply(W_match,length),sapply(W_match,length)) ),
                dims=c(N,N))
# Fixed parameters' values
rho <- 0.9
sigma <- 0.1
tau <- sigma^{-2}
# Realization of a zero-mean CAR distribution
Q <- tau*(D-rho*A)
S<-1
Q <- forceSymmetric(Q)
L <- Cholesky(Q,LDL=F)
z <- Matrix(rnorm(N*S),N,S)
u_per <- solve(L,z,system="Lt")
u <- solve(L,u_per,system="Pt")
# Realization of a covariate 
x1 <- runif(N,0,1)
# Linear combination for the log rate
eta <- 1 + 2*x1 + u
# Population at risk
Pop <- c(2454821,404589,301944,472298,368856,203686,246667,157973,1421781,343746,
         1767016,283548,308116,280813,287253,285642,242148,132572,278729,248480)
# Realization of a Poisson regression with rates Pop*exp(eta) 
mu_y <- numeric(N)
y <- numeric(N)
for(i in 1:N){
  mu_y[i] = Pop[i]*exp(eta[i])
  y[i] = rpois(1,mu_y[i])
}
# Maps of eta and u
counties$eta <- as.vector(eta)
plot(counties["eta"])
counties$u <- as.vector(u)
plot(counties["u"])

### Bayesian estimation
## Stan models - use the same priors in order to compare models
# Stan model using an explicit precision matrix for the CAR precision matrix 
car_full <- cmdstan_model("/home/x_romag/Documents/Tests/CAR_full.stan")
# Stan model using a sparse representation for the CAR precision matrix
car_sparse <- cmdstan_model("/home/x_romag/Documents/Tests/CAR_sparse.stan")
## Data inputs for the Stan models
# Some useful objects
X = as.matrix(x1)
P = ncol(X)
log_pop = log(Pop)
N_sparse = length(A@x)
A_w<-A@x
A_v<-A@i +1
A_u<-A@p +1
D_ii<-sapply(W_match,length)
log_det_D<-sum(log(D_ii))
lambda<-Schur(C, vectors = FALSE)$EValues
# Explicit model
data_full <- list(N = N,
                  P= P,
                  y = y,
                  log_offset = log_pop,
                  X = X,
                  A = as.matrix(A),
                  D = as.matrix(D),
                  lambda = lambda)
# Sparse model
data_sparse <- list(N = N,
                    P= P,
                    y = y,
                    log_offset = log_pop,
                    X = X,
                    N_sparse = N_sparse,
                    A_w = A_w,
                    A_v = A_v,
                    A_u = A_u,
                    D_ii = D_ii,
                    log_det_D = log_det_D,
                    lambda = lambda)
## MCMC sampling
# Explicit model
car_full_fit <- car_full$sample(data = data_full,
                                chains = 4,
                                refresh = 1000,
                                iter_warmup = 5000,
                                iter_sampling = 10000,
                                init = list(
                                  list(alpha = 1, beta = 2, phi_rho = 0.9, phi_tau =100),
                                  list(alpha = 1, beta = 2, phi_rho = 0.9, phi_tau =100),
                                  list(alpha = 1, beta = 2, phi_rho = 0.9, phi_tau =100),
                                  list(alpha = 1, beta = 2, phi_rho = 0.9, phi_tau =100)
                                )
)
# Sparse model
car_sparse_fit <- car_sparse$sample(data = data_sparse,
                                    chains = 4,
                                    refresh = 1000,
                                    iter_warmup = 5000,
                                    iter_sampling = 10000,
                                    init = list(
                                      list(alpha = 1, beta = 2, phi_rho = 0.9, phi_tau =100),
                                      list(alpha = 1, beta = 2, phi_rho = 0.9, phi_tau =100),
                                      list(alpha = 1, beta = 2, phi_rho = 0.9, phi_tau =100),
                                      list(alpha = 1, beta = 2, phi_rho = 0.9, phi_tau =100)
                                    )
)
## Posterior summaries
# Explicit model
car_full_fit$summary(variables = c("alpha", "beta", "phi_rho", "phi_sigma"),
                     posterior::default_summary_measures()
)
# Sparse model
car_sparse_fit$summary(variables = c("alpha", "beta", "phi_rho", "phi_sigma"),
                       posterior::default_summary_measures()
)
