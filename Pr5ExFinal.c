#include "FuncionsQR.h"

int Pr5ExFinal()
{
    // Configuración inicial para el RTBP
    int m = 3; // Dimensión del espacio de configuración (3 posiciones y 3 momentos)
    double Mu = 3.040357143e-6;
    double hh = -1.500384; // Nivel de energía

    // Condiciones iniciales para T y x
    double xx[7] = {3.051858, -0.988950, 0.0, 0.003235, 0.0, -0.999225, 0.0};
    double cg[7] = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Hiperplano definido por q1 = 0

    // Tolerancias
    double tol = 1e-10; // Para el método de Newton
    int maxit = 50; // Máximo número de iteraciones

    double h0rk = 0.1 , hminrk = 1e-6, hmaxrk = 0.5, tolrk = 1e-13; //Parámetros para flux()
    int npasmxrk = 1e5;

    Prm *Pmu=malloc(sizeof(Prm));
    Pmu->mu=Mu;

    // Llamada a opham para encontrar la órbita periódica
    int result = opham(m, hh, xx, tol, maxit, cg, h0rk, hminrk, hmaxrk, tolrk, npasmxrk, rtbp_h, rtbp, Pmu);

    if (result >= 0) {
        printf("\n\n\nPr5ExFinal: Conversion en %d iteraciones\n", result+1);
        printf("Orbita periodica encontrada:\n");
        printf("T = %.14f\n", xx[0]);
        printf("x = [%.16f, %.15f, %.18f, %.15f, %.16f, %.15f]\n", xx[1], xx[2], xx[3], xx[4], xx[5], xx[6]);
    } else {
        printf("Pr5ExFinal: No conversion en el numero maximo de iterados\n");
        printf("Last state:\n");
        printf("T = %.15f\n", xx[0]);
        printf("x = [%.15f, %.15f, %.15f, %.15f, %.15f, %.15f]\n", xx[1], xx[2], xx[3], xx[4], xx[5], xx[6]);
    }


    double t = 0.0, h = 0.1;
    double x_orb[6] = {xx[1], xx[2], xx[3], xx[4], xx[5], xx[6]};


    flux(&t, x_orb, &h, xx[0], hminrk, hmaxrk, tolrk, npasmxrk, 6, rtbp, Pmu, 55);

    printf("Orbit saved to orbHalo.txt.\n");

return result;
}

