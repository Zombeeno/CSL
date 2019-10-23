/*このプログラム:"show_complex"関数の導入で行列表示のコードが簡略化されたproto1.c*/
#include<math.h>
#include<stdio.h>
#include<complex.h>
#include<stdlib.h>
#include "cmatrix.h"
/*複素数はdouble complex型で定義でき、虚数単位はIで入力*/

void show_complex(double complex z){
  if(cimag(z)>=0){
  printf("%f+%fi ",creal(z),cimag(z));
  }
  else{
    printf("%f%fi ",creal(z),cimag(z));
  }
}

int main(void)
{
  double complex **A;
  A = alloc_cmatrix(2,2);
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      mat_elem(A,i,j)=i*1.0+j*1.0*I;
    }
  }
  printf("A=\n");
  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
    show_complex(mat_elem(A,i,j));}
    printf("\n");
  }
  
  return 0;
}
