#include "FuncionsQR.h"

int Pr5Ex4()
{
    int m = 1;
    double hh = -cos(0.5);
    double xx[3] = {2*M_PI, 0.5, 0};
    double cg[3] = {0.0, 1.0, 0.0};

    double tol = 1e-10;
    int maxit = 10;

    double hminrk = 1e-8, hmaxrk = 1e-4, tolrk = 1e-13;
    double h0rk = 1e-4;
    int npasmxrk = 1e6;

    Prm *P=malloc(sizeof(Prm));

    int iters = opham(m, hh, xx, tol, maxit, cg, h0rk, hminrk, hmaxrk, tolrk, npasmxrk, Ham, pendolvariacional,P);

    if (iters >= 0)
    {
        printf("\n\nPr5Ex4: Conversion en %d iteraciones\n", iters);
        printf("Orbita periodica encontrada:\n");
        printf("T = %f, q = %f, p = %f\n", xx[0], xx[1], xx[2]);
    }
    else
    {
        printf("Pr5Ex4: No conversion en el numero maximo de iteraciones %d\n", iters);
        printf("Last state:\n");
        printf("T = %f, q = %f, p = %f\n", xx[0], xx[1], xx[2]);
    }

return iters;
}
