#include "FuncionsQR.h"

int Pr1Ex2(void){
    int files=3, cols=2;
    double tol=1e-10, tk=0, sk=0;
    double **a, *x, *b, **r;

    printf("EXERCICI 2 PRACTICA 1\n");

    matriu_dinamica(&a, files, cols);
    vector_dinamic(&x, cols);
    vector_dinamic(&b, files);

    //a={{0,-4},{0,0},{5,-2}};
    a[0][0]=0;
    a[0][1]=-4;

    a[1][0]=0;
    a[1][1]=0;

    a[2][0]=5;
    a[2][1]=-2;

    b[0]=1;
    b[1]=3;
    b[2]=2;

    //b={1,3,2};

    QRsolve (files, cols+1, a, x, tol, b, &tk, &sk);

    r=amplia_matriu(files, cols, a, b);

    if (sk>=tol || sk<=-tol){
        system_solve(cols+1,r,x);
    }
    else return -1;

    esborra_matriu(a, files);
    free(x);
    free(b);
return 0;
}
