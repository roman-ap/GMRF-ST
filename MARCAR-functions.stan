/*
* Log probability density function of the Multivariate Autoregressive model
* of order 1 with CAR covariance/precision matrix (MARCAR): 
*
* @param y Spatio-Temporal process to model
* @param mu Mean vector
* @param rho_time Temporal dependence parameter
* @param rho_space Spatial dependence parameter
* @param sigma_space Scale parameter of the spatial dimension
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
* @return Log probability density
*/
real marcar_normal_lpdf ( vector y,
                          vector mu,
                          real rho_time,
                          real rho_space,
                          real sigma_space,
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
  
  vector[T*S] z = y - mu;
  real zItDsz;
  real zItAsz;
  real zAtDsz;
  real zAtAsz;
  real zIItDsz;
  real zIItAsz;
  vector[S] log_det_IsmrhoCs;
  zItDsz = dot_product(z .* ItDs_ii, z);
  zItAsz = dot_product(z, csr_matrix_times_vector(T*S, T*S, ItAs_w, ItAs_v, ItAs_u, z));
  zAtDsz = dot_product(z, csr_matrix_times_vector(T*S, T*S, AtDs_w, AtDs_v, AtDs_u, z));
  zAtAsz = dot_product(z, csr_matrix_times_vector(T*S, T*S, AtAs_w, AtAs_v, AtAs_u, z));
  zIItDsz = dot_product(z .* IItDs_ii, z);
  zIItAsz = dot_product(z, csr_matrix_times_vector(T*S, T*S, IItAs_w, IItAs_v, IItAs_u, z));
  for (s in 1:S){
    log_det_IsmrhoCs[s] = log1m( rho_space*lambdas[s] );
  }
  return 0.5 * (
    - T * S * log( 2 * pi() )
    + T * log1m( square(rho_time) )
    - 2 * square(S) * log( sigma_space )
    + S * log_det_Ds
    + S * sum ( log_det_IsmrhoCs )
    - inv_square(sigma_space) * ( zItDsz
                                  - rho_space*zItAsz
                                  - rho_time*zAtDsz
                                  + rho_time*rho_space*zAtAsz
                                  + square(rho_time)*zIItDsz
                                  - square(rho_time)*rho_space*zIItAsz ) 
    );
}

/*
* Log probability density function of the Multivariate Autoregressive model
* of order 1 with CAR covariance/precision matrix (MARCAR): 
*
* @param y Spatio-Temporal process to model
* @param mu Mean vector
* @param sigma_time Scale parameter of the temporal dimension
* @param sigma_space Scale parameter of the spatial dimension
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
* @return Log probability density
*/
/*
real marcar2_normal_lpdf ( vector y,
                           vector mu,
                           real sigma_time,
                           real sigma_space,
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
  
  vector[T*S] z = y - mu;
  real zItDsz;
  real zItAsz;
  real zAtDsz;
  real zAtAsz;
  real zIItDsz;
  real zIItAsz;
  vector[S] log_det_IsmrhoCs;
  zItDsz = dot_product(z .* ItDs_ii, z);
  zItAsz = dot_product(z, csr_matrix_times_vector(T*S, T*S, ItAs_w, ItAs_v, ItAs_u, z));
  zAtDsz = dot_product(z, csr_matrix_times_vector(T*S, T*S, AtDs_w, AtDs_v, AtDs_u, z));
  zAtAsz = dot_product(z, csr_matrix_times_vector(T*S, T*S, AtAs_w, AtAs_v, AtAs_u, z));
  zIItDsz = dot_product(z .* IItDs_ii, z);
  zIItAsz = dot_product(z, csr_matrix_times_vector(T*S, T*S, IItAs_w, IItAs_v, IItAs_u, z));
  for (s in 1:S){
    log_det_IsmrhoCs[s] = log1m( rho_space*lambdas[s] );
  }
  return 0.5 * (
    - T * S * log( 2 * pi() )
    - 2 * square(T) * log( sigma_time )
    + T * log1m( square(rho_time) )
    - 2 * square(S) * log( sigma_space )
    + S * log_det_Ds
    + S * sum ( log_det_IsmrhoCs )
    - inv_square(sigma_time * sigma_space) * ( zItDsz
                                               - rho_space*zItAsz
                                               - rho_time*zAtDsz
                                               + rho_time*rho_space*zAtAsz
                                               + square(rho_time)*zIItDsz
                                               - square(rho_time)*rho_space*zIItAsz ) 
    );
}
*/
