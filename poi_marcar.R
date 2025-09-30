#####  ##### 

### Set-up
library(sf)
library(dplyr)
library(Matrix)
library(cmdstanr)
library(posterior)
library(bayesplot)

### Sweden's regions
# Spatial geometries
precints <- read_sf("/proj/wildman/FixedKretsar/FixedKretsShape2023_2024.shp")
counties <- precints |> group_by(LANSKOD) |> summarise()
counties <- counties[-8,]
plot(counties)
# Number of counties
S<-nrow(counties)
# Binary spatial relationships  
S_match<-st_relate(counties,counties,pattern="****11212",sparse=TRUE)
# Sparse matrices for the CAR distribution
Ds_ii<-sapply(S_match,length)
D_space<-sparseMatrix(i=c( 1:S ),
                      j=c( 1:S ),
                      x=c( Ds_ii ),
                      dims=c(S,S))
A_space<-sparseMatrix(i=c( rep(1:S,Ds_ii) ),
                      j=c( unlist(S_match) ),
                      x=c( rep(1,length(unlist(S_match))) ),
                      dims=c(S,S))
C_space<-sparseMatrix(i=c( rep(1:S,Ds_ii) ),
                      j=c( unlist(S_match) ),
                      x=c( rep(1/Ds_ii,Ds_ii) ),
                      dims=c(S,S))
# Number of years
T<-21L
# Sparse matrices for the AR(1) distribution
I_time<-Diagonal(T, 1)
II_time<-sparseMatrix(i=c( 2:(T-1) ),
                      j=c( 2:(T-1) ),
                      x=c( rep(1, T-2) ),
                      dims=c(T,T))
A_time<-bandSparse(T, 
                   k = -1,
                   diagonals = list(rep(1, T-1)),
                   symmetric = TRUE)


# Fixed parameters' values
rho_time <- 0.7
rho_space <- 0.7
sigma_space <- 0.2
# Realization of a zero-mean CAR distribution
Q_time <- (I_time + rho_time^{2}*II_time - rho_time*A_time)
Q_space <- sigma_space^{-2}*(D_space - rho_space*A_space)
Q <- Matrix::kronecker(Q_time, Q_space)
M<-1
Q <- forceSymmetric(Q)
L <- Cholesky(Q,LDL=F)
z <- Matrix(rnorm(T*S*M),T*S,M)
u_per <- solve(L,z,system="Lt")
u <- solve(L,u_per,system="Pt")
# Linear combination for the log rate
eta <- u
# Population at risk
Pop <- rep(c(2454821,404589,301944,472298,368856,203686,246667,157973,1421781,343746,
             1767016,283548,308116,280813,287253,285642,242148,132572,278729,248480),
           T)
# Realization of a Poisson regression with rates Pop*exp(eta) 
mu_y <- numeric(T*S)
y <- numeric(T*S)
for(i in 1:(T*S)){
  mu_y[i] = Pop[i]*exp(eta[i])
  y[i] = rpois(1,mu_y[i])
}
inputs <- data.frame(y = y, log_pop = log(Pop))

### Bayesian estimation
# Stan model using a sparse representation for the CAR precision matrix
marcar_model <- cmdstan_model("/home/x_romag/Documents/GMRF-ST/Poisson_MARCAR.stan")
# Kronecker products
ItDs<-Matrix::kronecker(I_time, D_space)
ItAs<-Matrix::kronecker(I_time, A_space)
AtDs<-Matrix::kronecker(A_time, D_space)
AtAs<-Matrix::kronecker(A_time, A_space)
IItDs<-Matrix::kronecker(II_time, D_space)
IItAs<-Matrix::kronecker(II_time, A_space)
# Data inputs for the Stan models
##X = model.matrix(y ~ 1, data = inputs)
##P = ncol(X)
log_pop = inputs$log_pop
log_det_Ds<-sum(log(Ds_ii))
lambdas<-Schur(C_space, vectors = FALSE)$EValues
ItDs_ii = ItDs@x
N_ItAs = length(ItAs@x)
ItAs_w = ItAs@x
ItAs_v = ItAs@i + 1
ItAs_u = ItAs@p + 1
N_AtDs = length(AtDs@x)
AtDs_w = AtDs@x
AtDs_v = AtDs@i + 1
AtDs_u = AtDs@p + 1
N_AtAs = length(AtAs@x)
AtAs_w = AtAs@x
AtAs_v = AtAs@i + 1
AtAs_u = AtAs@p + 1
IItDs_ii = diag(IItDs)
N_IItAs = length(IItAs@x)
IItAs_w = IItAs@x
IItAs_v = IItAs@i + 1
IItAs_u = IItAs@p + 1
# Sparse model
marcar_data <- list(T = T,
                    S = S,
                    ##P= P,
                    y = y,
                    log_offset = log_pop,
                    ##X = X,
                    log_det_Ds = log_det_Ds,
                    lambdas = lambdas,
                    ItDs_ii = ItDs_ii,
                    N_ItAs = N_ItAs,
                    ItAs_w = ItAs_w,
                    ItAs_v = ItAs_v,
                    ItAs_u = ItAs_u,
                    N_AtDs = N_AtDs,
                    AtDs_w = AtDs_w,
                    AtDs_v = AtDs_v,
                    AtDs_u = AtDs_u,
                    N_AtAs = N_AtAs,
                    AtAs_w = AtAs_w,
                    AtAs_v = AtAs_v,
                    AtAs_u = AtAs_u,
                    IItDs_ii = IItDs_ii,
                    N_IItAs = N_IItAs,
                    IItAs_w = IItAs_w,
                    IItAs_v = IItAs_v,
                    IItAs_u = IItAs_u)
## MCMC sampling
# Sparse model
marcar_fit <- marcar_model$sample(data = marcar_data,
                                  chains = 4,
                                  parallel_chains = 4,
                                  refresh = 1000,
                                  iter_warmup = 5000,
                                  iter_sampling = 10000
                                  )
## Posterior summaries
# Sparse model
marcar_fit$summary(variables = c("rho_time", "rho_space", "sigma_space"))
marcar_fit$summary(variables = c("phi"), c("mean", "sd"), quantiles = ~ quantile2(., probs = c(0.05, 0.95)))

mcmc_hist(marcar_fit$draws("sigma_space"))
