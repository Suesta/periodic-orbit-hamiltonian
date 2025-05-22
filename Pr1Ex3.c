#include "FuncionsQR.h"


void Pr1Ex3(void){
    int files=3, cols=3;
    double tol=1e-10, tk=0, sk=0;
    double **a, *x, *b, **r;

    printf("\nEXERCICI 3 PRACTICA 1\n");

    matriu_dinamica(&a, files, cols);
    vector_dinamic(&x, cols);
    vector_dinamic(&b, files);

    a[0][0]=12;
    a[0][1]=-51;
    a[0][2]=4;

    a[1][0]=6;
    a[1][1]=167;
    a[1][2]=-68;

    a[2][0]=-4;
    a[2][1]=24;
    a[2][2]=-41;

    b[0]=1;
    b[1]=2;
    b[2]=3;

    QRsolve (files, cols+1, a, x, tol, b, &tk, &sk);

    r=amplia_matriu(files, cols, a, b);

    if (sk>=tol || sk<=-tol){
        system_solve(cols+1,r,x);
    }

    esborra_matriu(a, files);
    free(x);
    free(b);
}
