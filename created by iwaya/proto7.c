/*このプログラム:LAPACKのルーチンzgbsvを用いて、複素行列BについてのBX=Vなる解Xを与えるプログラム。Bは今回考える行列の対角成分をi倍したものをとり、Vは(0 0 0 34)とした(この時、結果が実部・虚部共に整数となる)*/
#include<math.h>
#include<stdio.h>
#include<complex.h>
#include<stdlib.h>
#include "cmatrix.h"
/*複素数はdouble complex型で定義でき、虚数単位はIで入力*/

/*zgbsv: A*X=B */
void zgbsv_(int *N, int *KL, int *KU, int *NRHS, double complex *AB, int *LDAB, int *IPIV, double complex *B, int *LDB, int *INFO);


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
    int N=5;
    int kl=2;
    int ku=2;
    int info;
    int *ipiv;
    int NRHS=1;
    int LDAB=2*kl+ku+1;
    double complex **A, **B, **C, **D, **AB;
    double complex *V;
    A = alloc_cmatrix(N,N);
    B = alloc_cmatrix(N,N);
    C = alloc_cmatrix(N,N);
    D = alloc_cmatrix(N,N);
    AB= alloc_cmatrix(LDAB,N);
    ipiv=alloc_ivector(N);
    V = alloc_cvector(N);
    for(int i=0;i<N;i++){
        vec_elem(V,i)=0;
        for(int j=0;j<N;j++){
            mat_elem(A,i,j) = 0;
            mat_elem(B,i,j) = 0;
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
    
    /*Dの生成*/
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            mat_elem(D,i,j)=3*mat_elem(A,i,j)+2*mat_elem(B,i,j)+mat_elem(C,i,j);
        }
    }
    
    /*ABの生成*/
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            if(mat_elem(D,i,j)!=0){
                mat_elem(AB,kl+ku+i-j,j)=mat_elem(D,i,j);
            }
        }
    }
    printf("matrix AB\n");
    print_cmatrix(LDAB,N,AB);
    
    
    /*Vの生成*/
    vec_elem(V,N-1)=34.0;
    printf("the original vector is\n");
    print_cvector(N,V);
    
    
    /*AX=Bの方程式の求解ルーチンzgbsvの引数の説明
     zgbsv(int N:Aの次元,int KL:Aの下帯の数,int KU:Aの上帯の数,int NRHS:Bのサイズ(N*NRHSのサイズ),double complex AB:行列Aの要素からなる「Aとは異なる」行列;出力時は変化,int LDAB:2*KL+KU+1以上の整数,int IPIV:ピボットを表す整数の配列, double complex B:入力時はB;出力時はX, int LDB:Nと同じで良い, int INFO:0でないならバグがある)
     AB,LDABについての詳細な記述はhttps://www.nag.co.uk/numeric/fl/nagdoc_fl26/pdf/f07/f07bnf.pdfを参照*/
    zgbsv_(&N, &kl, &ku, &NRHS, mat_ptr(AB), &LDAB, ipiv, V, &N, &info);
    /*ルーチンzgbsvにより、BX=Vの解XがVとして出力される。この時、ABも書き換えられてしまう*/
    printf("factorized AB\n");
    print_cmatrix(LDAB,N,AB);
    printf("the solution vector is\n");
    print_cvector(N,V);
    free_cmatrix(A);
    free_cmatrix(B);
    free_cmatrix(C);
    free_cmatrix(D);
    free_cmatrix(AB);
    free_cvector(V);
    return 0;
}
