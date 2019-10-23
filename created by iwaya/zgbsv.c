#include<math.h>
#include<stdio.h>
#include<complex.h>
#include<stdlib.h>
#include "cmatrix.h"
/*複素数はdouble complex型で定義でき、虚数単位はIで入力*/

int main(void)
{
  double complex **A;
  A = alloc_cmatrix(2,2);
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      mat_elem(A,i,j)=i*1.0+j*1.0*I;
    }
  }
        printf("A=\n%f+%f i %f+%f i\n%f+%f i %f+%f i \n",creal(mat_elem(A,0,0)),cimag(mat_elem(A,0,0)),creal(mat_elem(A,0,1)), cimag(mat_elem(A,0,1)),creal(mat_elem(A,1,0)), cimag(mat_elem(A,1,0)),creal(mat_elem(A,1,1)), cimag(mat_elem(A,1,1)));
        
  return 0;
}