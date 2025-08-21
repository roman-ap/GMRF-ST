/*
* Log probability density function of the conditional autoregressive (CAR) model: 
* WCAR specifications only.
*
* @param y Process to model
* @param mu Mean vector
* @param sigma Scale parameter
* @param rho Spatial dependence parameter
* @param A_w Sparse representation of the symmetric connectivity matrix, A
* @param A_v Column indices for values in A_w
* @param A_u Row starting indices for values in A_u
* @param D_ii The row sums of A
* @param log_det_D Log determinant of D = diag(D_ii)
* @param lambda Eigenvalues of C; for the WCAR specification, C is row-standardized.
* @param n Length of y
*
* @return Log probability density
*/
real wcar_normal_lpdf ( vector y,
                        vector mu,
                        real sigma,
                        real rho,
                        vector A_w, array[] int A_v, array[] int A_u,
                        vector D_ii, 
                        real log_det_D,
                        vector lambda,
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
    - 2 * n * log( sigma )
    + log_det_D
    + sum ( log_det_ImrhoC )
    - inv_square(sigma) * ( ztDz - rho*ztAz ) 
    );
}
