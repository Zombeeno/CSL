#include "matrix_util.h"
#include <stdio.h>
#include <math.h>

int main() {
    int n=100;
    int l;
    int m=100000;
    int i, j, k;
    double cx,cy;
    double **a;
    double **b;
    double **c;
    double **u;
    double **x;
    l = (n+1)*(n+1);
    printf("%d\n",b);

    /*初期値の設定*/
    a = alloc_dmatrix(n+1, n+1);
    b = alloc_dmatrix(n+1, n+1);
    c = alloc_dmatrix(n+1, n+1);
    for(i=0; i < n+1; ++i){
        for(j=0; j < n+1; ++j){
            if(i == j){
            a[i][j] = 1;
            }
            else{
            a[i][j]=0;
            }
        }
    }
    
    for(i=0; i < n; ++i){
        for(j=0; j < n+1; ++j){
            if(i == j){
                b[i][j] = -2;
            }
            else if(i == j+1){
                b[i][j] = 1;
            }
            else if(I == j-1){
                b[i][j] = 1;
            }
            b[i][j] = 0;
        }
    }
    for(j=0; j < n+1; ++j){
        if(i == j){
            b[n+1][j] = -2;
        }
        else if(i == j+1){
            b[n+1][j] = 2;
        }
         b[n+1][j] = 0;
        }
    
    
    fprint_dmatrix(stdout, n+1, n+1, a);
    
    /*Jacobi法の開始*/
    u = alloc_dmatrix(n+1,n+1);
    for(i = 0; i < m+1; ++i){
        for (j = 1; j < n; ++j){
            for (k = 1; k < n; ++k) {
                u[j][k] = (a[j+1][k]+a[j-1][k]+a[j][k+1]+a[j][k-1])/4;
                                     }
        }
            for (j = 1; j < n; ++j) {
                for (k = 1; k < n; ++k) {
                    a[j][k] = u[j][k];
                }
            }
    }
    fprint_dmatrix(stdout, n+1, n+1, u);
    /*gnuplot用*/
    x = alloc_dmatrix(b, 3);
    double n1;
    n1 = n;
    for(i = 0; i < n+1; ++i){
        for(j = 0; j < n+1; ++j){
            x[(n+1)*i+j][0] = i/n1;
            x[(n+1)*i+j][1] = j/n1;
            x[(n+1)*i+j][2] = a[i][j];
        }
    }
    fprint_dmatrix(stdout, b, 3, x);
    
    FILE *fp;
    char *filename = "Jacobi_solution.dat";
    fp = fopen(filename, "w");
    fprint_dmatrix(fp, b, 3, x);
    fclose(fp);
    free_dmatrix(a);
    free_dmatrix(u);
    free_dmatrix(x);
}
