#include "FuncionsQR.h"

int Pr4Ex2()
{
    double t = 0, h = 0.01, T = 100;
    double *x;
    double hmin=1e-4, hmax=0.01, tol=1e-8;
    int n=3, npasmax = 50000;

    Prm *PLorenz = malloc(sizeof(Prm));
    PLorenz->sigma = 10;
    PLorenz->rho = 28;
    PLorenz->beta = 8/3;

    x = malloc(n * sizeof(double));
    x[0] = 1;
    x[1] = 1;
    x[2] = 1;

    flux(&t, x, &h, T, hmin, hmax, tol, npasmax, n, lorenz_system, PLorenz,42);

    free(PLorenz);

    printf("Resultados guardados en lorenz_data1.dat\n");


    return 0;
}

