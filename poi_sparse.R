library(sf)
library(dplyr)
library(Matrix)
library(cmdstanr)
library(posterior)
library(bayesplot)


precints <- read_sf("/proj/wildman/FixedKretsar/FixedKretsShape2023_2024.shp")
counties <- precints |> group_by(LANSKOD) |> summarise()
counties <- counties[-8,]
plot(st_geometry(counties))

N<-nrow(counties)
W_match<-st_relate(counties,counties,pattern="****11212",sparse=TRUE)
D_ii<-sapply(W_match,length)
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

rho <- 0.9
sigma <- 1
tau <- sigma^{-2}

Q <- tau*(D-rho*A)
S<-1
Q <- forceSymmetric(Q)
L <- Cholesky(Q,LDL=F)
z <- Matrix(rnorm(N*S),N,S)
u_per <- solve(L,z,system="Lt")
u <- solve(L,u_per,system="Pt")

eta <- u

Pop <- c(2454821,404589,301944,472298,368856,203686,246667,157973,1421781,343746,
         1767016,283548,308116,280813,287253,285642,242148,132572,278729,248480)

y <- numeric(N)
mu_y <- numeric(N)
for(i in 1:N){
  mu_y[i] = Pop[i]*exp(eta[i])
  y[i] = rpois(1,mu_y[i])
}

A_w<-A@x
A_v<-A@i +1
A_u<-A@p +1
log_det_D<-sum(log(D_ii))
lambda<-Schur(C, vectors = FALSE)$EValues


counties$u <- as.vector(u)
plot(counties["u"])

car <- cmdstan_model("/home/x_romag/Documents/Stan/CAR_sparse.stan")

data_list <- list(N = N,
                  y = y,
                  log_offset = log(Pop),
                  D_ii = D_ii,
                  log_det_D = log_det_D,
                  lambda = lambda,
                  N_sparse = length(A_w),
                  A_w = A_w,
                  A_v = A_v,
                  A_u = A_u)

car_fit <- car$sample(data = data_list,
                      chains = 4,
                      refresh = 100,
                      iter_warmup = 5000,
                      iter_sampling = 10000
)

car_fit$summary(variables = c("phi_rho", "phi_sigma"),
                posterior::default_summary_measures()
)
