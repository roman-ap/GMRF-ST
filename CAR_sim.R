library(sf)
library(dplyr)
library(Matrix)

precints <- read_sf("/proj/wildman/FixedKretsar/FixedKretsShape2023_2024.shp")
plot(st_geometry(precints))


counties <- precints |> group_by(LANSKOD) |> summarise() |> st_cast("MULTIPOLYGON")
plot(st_geometry(counties))

plot(counties)

counties <- counties[-8,]
plot(st_geometry(counties))

# Dimensions
N<-nrow(counties)
# Precision matrix
W_match<-st_relate(counties,counties,pattern="F***1****",sparse=TRUE)
W<-sparseMatrix(i=c( 1:N,rep(1:N,sapply(W_match,length)) ),
                j=c( 1:N,unlist(W_match) ),
                x=c( sapply(W_match,length),rep(-1,length(unlist(W_match))) ),
                dims=c(N,N))

D<-sparseMatrix(i=c( 1:N ),
                j=c( 1:N ),
                x=c( sapply(W_match,length) ),
                dims=c(N,N))

A<-sparseMatrix(i=c( rep(1:N,sapply(W_match,length)) ),
                j=c( unlist(W_match) ),
                x=c( rep(1,length(unlist(W_match))) ),
                dims=c(N,N))

rho <- 0.9
sigma <- 1
tau <- sigma^{-2}

Q <- tau*(D-rho*A)

S<-100000
Q <- forceSymmetric(Q)
L <- Cholesky(Q,LDL=F)
z <- Matrix(rnorm(N*S),N,S)
u_per <- solve(L,z,system="Lt")
u <- solve(L,u_per,system="Pt")


means_u = rowMeans(u)
vars_u = apply(u,1,var)
Sigma = diag(solve(Q))
