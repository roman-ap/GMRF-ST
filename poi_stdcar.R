library(sf)
library(dplyr)
library(Matrix)
library(cmdstanr)
library(posterior)
library(bayesplot)


precints <- read_sf("/proj/wildman/FixedKretsar/FixedKretsShape2023_2024.shp")
counties <- precints |> group_by(LANSKOD) |> summarise()
counties <- counties[-8,]
plot(counties)

S<-nrow(counties)
S_match<-st_relate(counties,counties,pattern="****11212",sparse=TRUE)
Ds_ii<-sapply(S_match,length)
D<-sparseMatrix(i=c( 1:S ),
                j=c( 1:S ),
                x=c( Ds_ii ),
                dims=c(S,S))
A<-sparseMatrix(i=c( rep(1:S, Ds_ii) ),
                j=c( unlist(S_match) ),
                x=c( rep(1, length(unlist(S_match))) ),
                dims=c(S,S))
C<-sparseMatrix(i=c( rep(1:S, Ds_ii) ),
                j=c( unlist(S_match) ),
                x=c( rep(1/Ds_ii, Ds_ii) ),
                dims=c(S,S))

rho <- 0.9

Q <- (D-rho*A)
M<-1
Q <- forceSymmetric(Q)
L <- Cholesky(Q,LDL=F)
z <- Matrix(rnorm(S*M),S,M)
u_per <- solve(L,z,system="Lt")
u <- solve(L,u_per,system="Pt")

eta <- u

Pop <- c(2454821,404589,301944,472298,368856,203686,246667,157973,1421781,343746,
         1767016,283548,308116,280813,287253,285642,242148,132572,278729,248480)

y <- numeric(S)
mu_y <- numeric(S)
for(s in 1:S){
  mu_y[s] = Pop[s]*exp(eta[s])
  y[s] = rpois(1, mu_y[s])
}

A_w<-A@x
A_v<-A@i +1
A_u<-A@p +1
log_det_D<-sum(log(Ds_ii))
lambda<-Schur(C, vectors = FALSE)$EValues


car <- cmdstan_model("/home/x_romag/Documents/GMRF-ST/poi_stdcar.stan")

data_list <- list(S = S,
                  y = y,
                  log_offset = log(Pop),
                  D_ii = Ds_ii,
                  log_det_D = log_det_D,
                  lambda = lambda,
                  S_sparse = length(A_w),
                  A_w = A_w,
                  A_v = A_v,
                  A_u = A_u)

car_fit <- car$sample(data = data_list,
                      chains = 4,
                      parallel_chains = 4,
                      refresh = 100,
                      iter_warmup = 1000,
                      iter_sampling = 1000
)

car_fit$summary(variables = c("phi_rho"),
                posterior::default_summary_measures())
