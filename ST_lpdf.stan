/*
* Log probability density function of the standard conditional autoregressive model:
* WCAR specifications only.
*
* @param y Process to model
* @param mu Mean vector
* @param rho Spatial dependence parameter
* @param A_w Sparse representation of the symmetric adjacency matrix, A
* @param A_v Column indices for values in A_w
* @param A_u Row starting indices for values in A_u
* @param D_ii The row sums of A
* @param log_det_D Log determinant of D = diag(D_ii)
* @param lambda Eigenvalues of C; for the WCAR specification, C = D^{-1}*A
* @param n Length of y
*
* @return Log probability density function
*/
real std_wcar_lpdf ( vector y,
                     vector mu,
                     real rho,
                     real log_det_D,
                     vector lambda,
                     vector D_ii,
                     vector A_w, array[] int A_v, array[] int A_u,
                     int n ) {

  vector[n] z = y - mu;
  real ztDz;
  real ztAz;
  vector[n] log_det_ImrhoC;
  ztDz = dot_product(z .* D_ii, z);
  ztAz = dot_product(z, csr_matrix_times_vector(n, n, A_w, A_v, A_u, z));
  for (i in 1:n){
    log_det_ImrhoC[i] = log1m( rho*lambda[i] );
  }
  return 0.5 * (
    - n * log( 2 * pi() )
    + log_det_D
    + sum ( log_det_ImrhoC )
    - ( ztDz - rho*ztAz )
    );
}

/*
* Log probability density function of the standard Multivariate Autoregressive
* with WCAR covariance/precision matrix (MARCAR) model:
*
* @param y Spatio-Temporal process to model
* @param mu Mean vector
* @param rho_time Temporal dependence parameter
* @param rho_space Spatial dependence parameter
* @param log_det_Ds Log determinant of Ds = diag(Ds_ii)
* @param lambdas Eigenvalues of Cs = Ds^{-1}As
* @param ItDs_ii Diagonal elements of ItDs = It⊗Ds
* @param ItAs_w Sparse representation of ItAs = It⊗As
* @param ItAs_v Column indices for values in ItAs
* @param ItAs_u Row starting indices for values in ItAs
* @param AtDs_w Sparse representation of AtDs = At⊗Ds
* @param AtDs_v Column indices for values in AtDs
* @param AtDs_u Row starting indices for values in AtDs
* @param AtAs_w Sparse representation of AtAs = At⊗As
* @param AtAs_v Column indices for values in AtAs
* @param AtAs_u Row starting indices for values in AtAs
* @param IItDs_ii Diagonal elements of IItDs = IIt⊗Ds
* @param IItAs_w Sparse representation of IItAs = IIt⊗As
* @param IItAs_v Column indices for values in IItAs
* @param IItAs_u Row starting indices for values in IItAs
* @param T Number of temporal units
* @param S Number of spatial units
*
* @return Log probability density function
*
*/
real std_marcar_lpdf ( vector y,
                       vector mu,
                       real rho_time,
                       real rho_space,
                       real log_det_Ds,
                       vector lambdas,
                       vector ItDs_ii,
                       vector ItAs_w, array[] int ItAs_v, array[] int ItAs_u,
                       vector AtDs_w, array[] int AtDs_v, array[] int AtDs_u,
                       vector AtAs_w, array[] int AtAs_v, array[] int AtAs_u,
                       vector IItDs_ii,
                       vector IItAs_w, array[] int IItAs_v, array[] int IItAs_u,
                       int T,
                       int S) {

  vector[T*S] u = y - mu;
  real utItDsu;
  real utItAsu;
  real utAtDsu;
  real utAtAsu;
  real utIItDsu;
  real utIItAsu;
  vector[S] log_det_IsmrhoCs;
  utItDsu = dot_product(u .* ItDs_ii, u); // t(u) %*% (I_T %x% D_S) %*% u
  utItAsu = dot_product(u, csr_matrix_times_vector(T*S, T*S, ItAs_w, ItAs_v, ItAs_u, u)); // t(u) %*% (I_T %x% A_S) %*% u
  utAtDsu = dot_product(u, csr_matrix_times_vector(T*S, T*S, AtDs_w, AtDs_v, AtDs_u, u)); // t(u) %*% (A_T %x% D_S) %*% u
  utAtAsu = dot_product(u, csr_matrix_times_vector(T*S, T*S, AtAs_w, AtAs_v, AtAs_u, u)); // t(u) %*% (A_T %x% D_S) %*% u
  utIItDsu = dot_product(u .* IItDs_ii, u); // t(u) %*% (I_T^{*} %x% D_S) %*% u
  utIItAsu = dot_product(u, csr_matrix_times_vector(T*S, T*S, IItAs_w, IItAs_v, IItAs_u, u)); // t(u) %*% (I_T_{*} %x% A_S) %*% u
  for (s in 1:S){
    log_det_IsmrhoCs[s] = log1m( rho_space*lambdas[s] );
  }
  return 0.5 * (
    - T * S * log( 2 * pi() )
    + T * log1m( square(rho_time) )
    + S * log_det_Ds
    + S * sum ( log_det_IsmrhoCs )
    - ( utItDsu - rho_space*utItAsu
        + square(rho_time)*utIItDsu - square(rho_time)*rho_space*utIItAsu
        - rho_time*utAtDsu + rho_time*rho_space*utAtAsu )
    );
}
