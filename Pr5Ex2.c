#include "FuncionsQR.h"

int Pr5Ex2() {

    int m = 1;
    int n=2*m;
    double hh = -cos(0.5);
    double xx[3] = {0.1, 0.5, 0.0};
    double cg[3] = {0.0, 1.0, 0.0};
    double hminrk = 1e-8, hmaxrk = 1e-3, tolrk = 1e-13;
    double h0rk = hmaxrk;
    int npasmxrk = 1e6;
    double *ff , **df;

    Prm *P=malloc(sizeof(Prm));
    vector_dinamic(&ff,n+2);
    matriu_dinamica(&df, n+2, n+1);

    int result = opham_fdf(m, hh, xx, cg, ff, df, h0rk, hminrk, hmaxrk, tolrk, npasmxrk, Ham, pendolvariacional, P);

    if (result != 0) {
        fprintf(stderr, "Error en opham_fdf\n");
        free(ff);
        esborra_matriu(df,n+2);
        return -1;
    }

    // Print results
    printf("Resultados de opham_fdf:");
    printf("\n\nF(X):\n");
    for (int i = 0; i < n+2; i++) {
        printf("%f\n", ff[i]);
    }

    printf("\n\nDF(X):\n\n");
    for (int i = 0; i < n+2; i++) {
        for (int j = 0; j < n+1; j++) {
            printf("%f\t ", df[i][j]);
        }
        printf("\n\n");
    }

    free(ff);
    esborra_matriu(df,n+2);

    return 0;
}

