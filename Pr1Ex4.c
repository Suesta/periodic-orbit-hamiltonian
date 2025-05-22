#include "FuncionsQR.h"


void Pr1Ex4(void){

    int n=0, contador=0;
    double tol=1e-10, tk=0, sk=0;
    double **a, **mx, *x, *b, **r;
    double MaxError=0;

    printf("\nEXERCICI 4 PRACTICA 1\n");

    printf("\nINTRODUEIX DIMENSIO: n=");
    scanf("%d", &n);

    matriu_dinamica(&a, n,n);
    matriu_dinamica(&mx, n,1);
    vector_dinamic(&x, n);
    vector_dinamic(&b, n);

    printf("\n\n\nMATRIU MX (VECTOR X DE 1'S)\n\n");
    for(int j=0;j<n;j++){
        mx[j][0]=1;
        printf("%g\n", mx[j][0]);
    }

    printf("\n\n\nMATRIU A\n\n");
    srand(time(0));
    for(int i=0;i<n;i++){
            printf("\n\n");
        for(int j=0;j<n;j++){
            a[i][j]=((double)rand())/RAND_MAX;
            printf("%g\t",a[i][j]);
        }
    }

    printf("\n\n\nVECTOR B ADIENT\n\n");
    for(int j=0;j<n;j++){
        b[j]=prod_matrius(n,n,1,a,mx)[j][0];
        printf("%g\t",b[j]);
    }
    printf("\n\n\n\n\n");

    QRsolve (n, n+1, a, x, tol, b, &tk, &sk);

    r=amplia_matriu(n, n, a, b);

    if (sk>=tol || sk<=-tol){
        system_solve(n+1,r,x);
    }
    while(contador<n){
        if(abs(x[contador]-1)>MaxError){
            MaxError=x[contador]-1;
            contador++;
        }
        else contador++;
    }

    printf("\n\nMaxError=%g\n\n", MaxError);

    esborra_matriu(a, n);
    free(x);
    free(b);
}
