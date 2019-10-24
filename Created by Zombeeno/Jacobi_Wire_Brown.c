#include "matrix_util.h"
#include <stdio.h>
#include <math.h>

double Jacobi(int l, double **D)
{
    int n=10000, i1, j1, k1;
    /*Jacobi法に必要なツール*/
    double error=pow(10,-5.0);
    double *sum;
    double *storage;
    double sum2;
    double sol;
    double *f;
    double *x;
    f = alloc_dvector(l+1);
    x = alloc_dvector(l+1);
    sum = alloc_dvector(l+1);
    storage = alloc_dvector(l+1);
    /*外力f*/
    for(i1=0; i1 < l+1; ++i1){
        if(i1 == l){
            f[i1] = 1;
        }
        f[i1] = 0;
       }
    for(k1=0; k1 < n; ++k1){
        for(i1=0; i1 < l+1; ++i1){
            storage[i1]=x[i1];
        }
        for(i1=0; i1 < l+1; ++i1){
            for(j1=0; j1 < l+1; ++j1){
                sum[i1] = sum[i1] + D[i1][j1]*x[j1];
                }
        }
        for(i1=0; i1 < l+1; ++i1){
            x[i1]=(f[i1]-sum[i1])/D[i1][i1];
            sum[i1] = 0;
        }
        for(i1=0; i1 < l+1; ++i1){
            sum2 = sum2 + x[i1]-storage[i1];
        }
        if(sum2 < error){
            break;
        }
    }
    sol = x[l];
free_dvector(f);
free_dvector(sum);
free_dvector(storage);
free_dvector(x);
    return sol;
}



int main(void) {
    int l;
    double dz;
    double d; /*mm*/
    double rho=19.25;/*g/cm^3*/
    double ReE;
    double ImE;
    double I;
    double L=600; /*mm*/
    double **w;
    double S=0.003311;
    int i, j, k;
    double **a;
    double **b;
    double **c;
    double **ReD;
    double **ImD;

    //要素数と初期値の代入および計算
    scanf("%d",&l);
    scanf("%lf",&d);
    dz=L/(2*l);
    I=rho/1000*L*pow(d,4.0)/2;
    printf("%d %lf %lf %lf",l,d,I,dz);
    
    //初期値の設定
    a = alloc_dmatrix(l+1, l+1);
    b = alloc_dmatrix(l+1, l+1);
    c = alloc_dmatrix(l+1, l+1);
    ReD = alloc_dmatrix(l+1, l+1);
    ImD = alloc_dmatrix(l+1, l+1);
    w = alloc_dmatrix(l+1, 2);
    
    //単位行列A
    for(i=0; i < l+1; ++i){
        for(j=0; j < l+1; ++j){
            if(i == j){
            a[i][j] = 1;
            }
            else{
            a[i][j]=0;
            }
        }
    }
    //行列B
    for(i=0; i < l; ++i){
        for(j=0; j < l+1; ++j){
            if(i == j){
                b[i][j] = -2;
            }
            else if(i == j+1){
                b[i][j] = 1;
            }
            else if(i == j-1){
                b[i][j] = 1;
            }
            b[i][j] = 0;
        }
    }
    for(j=0; j < l+1; ++j){
        if(j == l){
            b[l][j] = -2;
        }
        else if(j == l-1){
            b[l][j] = 2;
        }
         b[l][j] = 0;
        }
    //行列C
    for(i=0; i < l-1; ++i){
        for(j=0; j < l+1; ++j){
            if(i == j){
                c[i][j] = 6;
            }
            else if(i == j+1){
                c[i][j] = -4;
            }
            else if(i == j-1){
                c[i][j] = -4;
            }
            else if(i == j+1){
                c[i][j] = 1;
            }
            else if(i == j-1){
                c[i][j] = 1;
            }
            c[i][j] = 0;
        }
    }
    for(j=0; j < l+1; ++j){
        if(j == l-3){
            c[l-1][j] = 1;
        }
        else if(j == l-2){
            c[l-1][j] = -4;
        }
        else if(j == l-1){
            c[l-1][j] = 7;
        }
        else if(j == l){
            c[l-1][j] = -4;
        }
        c[l-1][j] = 0;
    }
    for(j=0; j < l+1; ++j){
        if(j == l-2){
            c[l][j] = 2;
        }
        else if(j == l-1){
            c[l][j] = -8;
        }
        else if(j == l){
            c[l][j] = 6;
        }
        c[l-1][j] = 0;
    }

    //行列Dの計算
    for(k=10; k < 10000; ++k){
    w[0][k] = k;
    for(i=0; i < l-1; ++i){
        for(j=0; j < l+1; ++j){
            ReD[i][j]=-ReE*I/pow(dz,4.0)*c[i][j]+S/pow(dz,2)*b[i][j]+rho*pow(w[0][k],2.0)*a[i][j];
            ImD[i][j]=-ImE*I/pow(dz,4.0)*c[i][j]+S/pow(dz,2)*b[i][j]+rho*pow(w[0][k],2.0)*a[i][j];
     }
    }
    w[1][k] = Jacobi(l,ReD);
    w[2][k] = Jacobi(l,ImD);
}
//ファイルへの書き込み
    FILE *fp;
    char *filename = "Wire_Brown.dat";
    fp = fopen(filename, "w");
    fprint_dmatrix(fp, l, 2, w);
    fclose(fp);
//容量の解放
    free_dmatrix(a);
    free_dmatrix(b);
    free_dmatrix(c);
    free_dmatrix(w);
    free_dmatrix(ReD);
    free_dmatrix(ImD);

return 0;

}
