/*このプログラム:複素行列A,B,Cを複素型で生成する。次数は34行目により可変*/
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

static inline void print_real(int m, int n, double complex **mat) {
  int i, j;
  printf("column:%d row:%d\n", m, n);
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j){printf("%f ",creal(mat_elem(mat,i,j)));}
    printf("\n");
  }
}

int main(void)
{
  /*N:生成する行列の次数、変更可能*/
  int N=5;
  double complex **A, **B, **C;
  A = alloc_cmatrix(N,N);
  B = alloc_cmatrix(N,N);
  C = alloc_cmatrix(N,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      mat_elem(A,i,j) = 0;
      mat_elem(B,i,j) = 0;
      mat_elem(C,i,j) = 0;
    }
  }
  
  /*Aの生成*/
  for(int i=0;i<N;i++){
    mat_elem(A,i,i)=1.0;
  }
  
  /*Bの生成*/
  mat_elem(B,0,0)=-2.0;
  mat_elem(B,0,1)= 1.0;
  for(int i=1;i<N-1;i++){
    mat_elem(B,i,i-1)=1.0;
    mat_elem(B,i,i)=-2.0;
    mat_elem(B,i,i+1)=1.0;
  }
  mat_elem(B,N-1,N-2)=2.0;
  mat_elem(B,N-1,N-1)=-2.0;
  
  /*Cの生成*/
  for(int i=0;i<N-2;i++){
    mat_elem(C,i,i+2) = 1.0;
    mat_elem(C,i+2,i) = 1.0;
  }
  for(int i=0;i<N-1;i++){
    mat_elem(C,i,i+1) = -4.0;
    mat_elem(C,i+1,i) = -4.0;
  }
  for(int i=0;i<N;i++){
    mat_elem(C,i,i) = 6.0;
  }
  mat_elem(C,N-2,N-2)= 7.0;
  mat_elem(C,N-1,N-3)= 2.0;
  mat_elem(C,N-1,N-2)=-8.0;
  
  /*A,B,Cの表示:現段階ではこれらは複素でないため、実部しか表示しない関数を使用*/
  printf("matrix A ");
  print_real(N,N,A);
  printf("matrix B ");
  print_real(N,N,B);
  printf("matrix C ");
  print_real(N,N,C);
  
  free_cmatrix(A);
  free_cmatrix(B);
  free_cmatrix(C);
  return 0;
}
