#include "Matrix.h"

void qr_comp(const Matrix<double>& A, Matrix<double>& Q, Matrix<double>& R){
  /* Input: square matrix A
   * output: orthogonal matrix Q, upper triangular matrix R with non-negative diagonal entries 
   * such that A = QR. */
  double a, sin_t, cos_t, s;
  int n = A.nrows();
  assert (A.ncols() == n && Q.nrows() == n && Q.ncols() == n && R.nrows() == n && R.ncols() == n);

  R.copy(A);  // initialize R = A
  Q.fill(0);
  for(int i = 0;i < n;i++) Q(i,i) = 1.0;
  
  for(int i = 0;i < n;i++){
    for(int j = i+1; j <n; j++){ // zero the remainder of column i with Givens rotations
      if(R(j,i) == 0) continue;
      s = ( (R(j,i) < 0) ?  -1.0 : 1.0);  
      if(R(i,i) > s*R(j,i)){// begin numerical hygene for rotating [r_ii,r_ji] -> [r_ii/sin_t, 0]
        a = R(j,i)/R(i,i);
        sin_t = 1/sqrt(1+a*a);
        cos_t = a*sin_t;
        R(i,i) /= sin_t; // result of applying rotation to (R(i,i),R(j,i))
      }
      else {
        a = R(i,i)/R(j,i); // continue hygene (case 2)
        cos_t = s/sqrt(1+a*a);
        sin_t = a*cos_t;
        R(i,i) = R(j,i)/cos_t; // result of applying rotation to (R(i,i),R(j,i))
      }
      R(j,i) = 0;
      for(int k = i+1;k < n;k++) { // apply rotation to the remainder of the entries in rows i, j
        a = sin_t*R(i,k) + cos_t*R(j,k);
        R(j,k) = -cos_t*R(i,k) + sin_t*R(j,k);
        R(i,k) = a;
      }
      for(int k = 0;k < n;k++){ // rotate the Q-matrix
        a = sin_t*Q(i,k) + cos_t*Q(j,k);
        Q(j,k) = -cos_t*Q(i,k) + sin_t*Q(j,k);
        Q(i,k) = a;
      }
    }// end j-loop
  }// end i-loop
}
