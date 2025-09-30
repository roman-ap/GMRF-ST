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
plot(st_geometry(counties))

### Simulated data for Swedish counties
# Number of counties
N<-nrow(counties)
# Binary spatial relationships  
W_match<-st_relate(counties,counties,pattern="****11212",sparse=TRUE)
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
rho <- 0.5
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
lambda <- numeric(N)
y <- numeric(N)
for(i in 1:N){
  lambda[i] = Pop[i]*exp(eta[i])
  y[i] = rpois(1,lambda[i])
}

counties$y <- as.vector(y)
counties$x1 <- x1
counties$u <- as.vector(u)
counties$Pop <- Pop 
counties$lambda <- as.vector(lambda)
plot(counties["u"])

### Bayesian estimation
# Stan model using an explicit precision matrix for the CAR modeling 
car <- cmdstan_model("/home/x_romag/Desktop/CAR_full.stan")
# Data inputs for the Stan model
data_list <- list(N = N,
                  P= 1,
                  X = as.matrix(x1),
                  y = y,
                  log_offset = log(Pop),
                  A = as.matrix(A),
                  D = as.matrix(D))
# MCMC sampling
car_fit <- car$sample(data = data_list,
                      chains = 4,
                      refresh = 1000,
                      iter_warmup = 5000,
                      iter_sampling = 10000
)
# Posterior summaries
car_fit$summary(variables = c("alpha", "beta", "rho", "phi_tau"),
                posterior::default_summary_measures()
)
