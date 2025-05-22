#include "FuncionsQR.h"


int Pr4Ex3(void)
{
    //double hmin = 1e-8, hmax = 1e-4, tol = 1e-8, delta = 1e-4; // Parámetros de integración ajustados
    //double hmin = 1e-6, hmax = 0.1, tol = 1e-8; // Parámetros de integración ajustados
    double hmin = 1e-10, hmax = 1e-2, tol = 1e-8;
    int n = 2, npasmax = 1e6; // Dimensión del sistema
    double t = 0.0, h = 1e-6, T = 0.5; // Intervalo de integración
    double *x;
    double *f;

    f = malloc((n+n*n) * sizeof(double));

    Prm *P = malloc(sizeof(Prm));

    P->a = 0.5;

    x = malloc((n+n*n) * sizeof(double));

    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 1.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 1.0;

    variacionales4(n+n*n,t,x,f,P);

    flux(&t, x, &h, T, hmin, hmax, tol, npasmax, n+n*n, variacionales4, P,43);

    printf("Resultados guardados en pr4ex3_results.dat\n");

    free(P);
    return 0;
}



