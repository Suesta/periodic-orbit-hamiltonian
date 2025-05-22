#include "FuncionsQR.h"

void E2(void);

void E2(void){
    int m=3, n=3, maxit=50;
    double tol=1e-12;
    double *x, *f, **df;

    vector_dinamic(&x,n);
    vector_dinamic(&f,m);
    matriu_dinamica(&df, m, n);

    //Llavor inicial
    x[0]=1;
    x[1]=1;
    x[2]=3;

    newton(m, n, x, fdf, maxit, tol, f, df);

}

