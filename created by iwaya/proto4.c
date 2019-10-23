/*このプログラム:Jacobi法のプロトタイプとして、複素行列BについてのBX=Vなる解Xを与えるプログラム。Bは今回考える行列の対角成分をi倍したものをとり、
Vは(0 0 14)とした(この時、結果が整数となる)*/
#include<math.h>
#include<stdio.h>
#include<complex.h>
#include<stdlib.h>
#include "cmatrix.h"
/*複素数はdouble complex型で定義でき、虚数単位はIで入力*/

/*複素数を表示する関数*/
void show_complex(double complex z){
  if(cimag(z)>=0){
    printf("%f+%fi ",creal(z),cimag(z));
  }
  else{
    printf("%f%fi ",creal(z),cimag(z));
  }
}

static inline void print_cvector(int n, double complex *vec) {
  int i;
  printf("components:%d\n", n);
  for (i = 0; i < n; ++i) {
    show_complex(vec_elem(vec, i));
    printf("\n");
  }
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
  int N=3;
  double complex **B;
  double complex *X,*Y,*V;
  double complex sum = 0;
  double gap = 0;
  double gsum = 0;
  B = alloc_cmatrix(N,N);
  X = alloc_cvector(N);
  Y = alloc_cvector(N);
  V = alloc_cvector(N);
  for(int i=0;i<N;i++){
    vec_elem(X,i)=0;
    vec_elem(Y,i)=0;
    vec_elem(V,i)=0;
    for(int j=0;j<N;j++){
      mat_elem(B,i,j) = 0;
    }
  }
  
  /*Bの生成*/
  mat_elem(B,0,0)=-2.0*I;
  mat_elem(B,0,1)= 1.0;
  for(int i=1;i<N-1;i++){
    mat_elem(B,i,i-1)=1.0;
    mat_elem(B,i,i)=-2.0*I;
    mat_elem(B,i,i+1)=1.0;
  }
  mat_elem(B,N-1,N-2)=2.0;
  mat_elem(B,N-1,N-1)=-2.0*I;
  printf("the Matrix\n");
  print_cmatrix(N,N,B);
  /*Vの生成*/
  vec_elem(V,N-1)=14.0;
  printf("the original vector is\n");
  print_cvector(N,V);
  
  /*Jacobi法:初期値は0-vector*/
  for(int a=0;a<10000;a++){
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){sum = sum+mat_elem(B,i,j)*vec_elem(X,j);}
      vec_elem(Y,i)=(vec_elem(V,i)-sum+mat_elem(B,i,i)*vec_elem(X,i))/mat_elem(B,i,i);
      gap = cabs((vec_elem(X,i)-vec_elem(Y,i)));
      gsum = gsum + gap;
      sum = 0;
    }
    if(gsum<0.000000001){
      printf("iteration count: %d\n",a);
    break;}
    for(int i=0;i<N;i++){vec_elem(X,i)=vec_elem(Y,i);}
    sum = 0;
    gsum = 0;
  }
  
  printf("the solution vector is\n");
  print_cvector(N,Y);
  
  free_cmatrix(B);
  free_cvector(V);
  free_cvector(X);
  free_cvector(Y);
  return 0;
}
