#####  ##### 

### Set-up
library(sf)
library(dplyr)
library(Matrix)

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
D_space<-sparseMatrix(i=c( 1:S ),
                      j=c( 1:S ),
                      x=c( sapply(S_match,length) ),
                      dims=c(S,S))
A_space<-sparseMatrix(i=c( rep(1:S,sapply(S_match,length)) ),
                      j=c( unlist(S_match) ),
                      x=c( rep(1,length(unlist(S_match))) ),
                      dims=c(S,S))
C_space<-sparseMatrix(i=c( rep(1:S,sapply(S_match,length)) ),
                      j=c( unlist(S_match) ),
                      x=c( rep(1/sapply(S_match,length),sapply(S_match,length)) ),
                      dims=c(S,S))
# Number of years
T<-10L
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
# Kronecker products
ItDs<-Matrix::kronecker(I_time, D_space)
ItAs<-Matrix::kronecker(I_time, A_space)
AtDs<-Matrix::kronecker(A_time, D_space)
AtAs<-Matrix::kronecker(A_time, A_space)
IItDs<-Matrix::kronecker(II_time, D_space)
IItAs<-Matrix::kronecker(II_time, A_space)
  
  

rho_space <- 0.9
sigma_space <- 1
rho_time <- 0.9
sigma_time <- 1
Q_time <- sigma_time^{-2}*(I_time + rho_time^{2}*II_time - rho_time*A_time)
Q_space <- sigma_space^{-2}*(D_space - rho_space*A_space)

Q <- Matrix::kronecker(Q_time, Q_space)

M<-10000
Q <- forceSymmetric(Q)
L <- Cholesky(Q,LDL=F)
z <- Matrix(rnorm(T*S*M),T*S,M)
u_per <- solve(L,z,system="Lt")
u <- solve(L,u_per,system="Pt")


means_u = rowMeans(u)
vars_u = apply(u,1,var)
Sigma = diag(solve(Q))
