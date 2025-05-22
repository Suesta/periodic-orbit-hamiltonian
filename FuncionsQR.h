#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define RTBP_M 3
#define RTBP_N 6
#define RTBP_NV1 42

#define SQR(x) ((x)*(x))

#define N RTBP_N
#define NV1 RTBP_NV1

#define X x[0]
#define Y x[1]
#define Z x[2]
#define PX x[3]
#define PY x[4]
#define PZ x[5]

typedef struct parametres{
    double tMax;
    double a;
    double sigma;
    double rho;
    double beta;
    double mu;
} Prm;

int Pr1Ex2(void);
void Pr1Ex3(void);
void Pr1Ex4(void);

void E2(void);

int Pr3Ex2(void);

int Pr4Ex2(void);
int Pr4Ex3(void);

int Pr5Ex2(void);
int Pr5Ex4(void);
int Pr5ExFinal(void);

//practica 1
int QRsolve (int m, int n, double **matrix, double *x, double tol, double *b, double *tk, double *sk);
void matriu_dinamica(double ***m, int files, int cols);
void esborra_matriu(double **m, int files);
void vector_dinamic(double **vector, int size);
double prod_vectors (int size, double *vf, double *vc);
double **amplia_matriu(int f, int c, double **m, double *b);
double **prod_matrius (int f, int fc, int c, double **a, double **b);
void system_solve (int c, double **m, double *x);

//practica 2
int fdf(int m, int n, double *x, double *f, double **df);
int newton(int m, int n, double *x,int (*fdf)(int m, int n, double *x, double *f, double **df),int maxit, double tol, double *f, double **df);

//practica 3
void camp1(int n, double t, double *x, double *f, Prm *prm);
void solucioExacta(double t, double x0, double y0, double *xExacta, Prm *p);

int rk45(double *t, double *x,  double *h, double hmin, double hmax, double tol, int n, void (*camp)(int n, double t, double *x, double *f, Prm *prm), Prm *prm);
int rk78(double *t, double *x,  double *h, double hmin, double hmax, double tol, int n, void (*camp)(int n, double t, double *x, double *f, Prm *prm), Prm *prm);

//practica 4
void lorenz_system(int n, double t, double *x, double *f, Prm *prm);
void variacionales4(int n, double t, double *x, double *f, Prm *prm);
int flux(double *t, double *x, double *h, double T, double hmin, double hmax, double tol, int npasmax, int n, void (*camp)(int n, double t, double *x, double *f,Prm *prm), Prm *prm, int identif);

//practica 5
int Ham(double *x, double *H, Prm *prm);
void pendolvariacional(int n, double t, double *x, double *f, Prm *prm);

int opham_fdf(int m, double hh, double *xx, double *cg, double *ff, double **df, double h0rk, double hminrk, double hmaxrk, double tolrk, int npasmxrk,int (*ham)(double *x, double *H, Prm *prm),void (*camp)(int n, double t, double *x, double *f, Prm *prm),Prm *prm);

int newton1(int (*opham_fdf)(int m, double hh, double *xx, double *cg, double *ff, double **df, double h0rk, double hminrk, double hmaxrk, double tolrk, int npasmxrk, int (*ham)(double *x, double *H, Prm *prm), void (*camp)(int n, double t, double *x, double *f, Prm *prm), Prm *prm),int m, double hh, double *xx, double *cg, double *ff, double **df, double h0rk, double hminrk, double hmaxrk, double tolrk, int npasmxrk,int (*ham)(double *x, double *H, Prm *prm),void (*camp)(int n, double t, double *x, double *f, Prm *prm),int maxit, double tol, Prm *prm);

int opham(int m, double hh, double *xx, double tol, int maxit, double *cg,double h0rk, double hminrk, double hmaxrk, double tolrk, int npasmxrk,int (*ham)(double *x, double *H, Prm *prm),void (*camp)(int n, double t, double *x, double *f, Prm *prm),Prm *prm);

void rtbp(int n, double t, double *x, double *f, Prm *prm);
int rtbp_h (double *x, double *H, Prm *prm);
