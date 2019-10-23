/*このプログラム:複素行列をプリントする関数"print_cmatrix"を導入した。*/
#include<math.h>
#include<stdio.h>
#include<complex.h>
#include<stdlib.h>
#include "cmatrix.h"
/*複素数はdouble complex型で定義でき、虚数単位はIで入力*/

void show_complex(double complex z){
  printf("%f+%fi ",creal(z),cimag(z));
}

static inline void print_cmatrix(int m, int n, double complex **mat) {
  int i, j;
  printf("column:%d row:%d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) show_complex(mat_elem(mat, i, j));
    printf("\n");
  }
}

int main(void)
{
  double complex **A;
  A = alloc_cmatrix(3,2);
  for(int i=0;i<3;i++){
    for(int j=0;j<2;j++){
      mat_elem(A,i,j)=i*1.0+j*1.0*I;
    }
  }
  
  print_cmatrix(3,2,A);
  return 0;
}
