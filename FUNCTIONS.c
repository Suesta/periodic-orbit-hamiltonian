#include "FuncionsQR.h"

int QRsolve (int m, int n, double **matrix, double *x, double tol, double *b, double *tk, double *sk)
{

    if(b==NULL)return -1;
    else{
        if(m<n-1)return -1;
        else{
            int k=1, cols=0, loop=0;
            double *a1, *uk, **a;
            double pivot=0, tau=0;

            vector_dinamic(&a1, m);
            vector_dinamic(&uk, m);
            a=amplia_matriu(m,n-1,matrix,b);
                /*printf("\n\nmatriu A iteracio %d\n\n", cols+1);
                for(int i=0;i<m;i++){
                        printf("\n\n");
                    for(int j=0;j<n;j++){
                        printf("%.3f\t",a[i][j]);
                    }
                }*/
            if (m==n-1)loop=n-1;
            else loop = n;

            while(k<loop){
                (*sk)=0;
                for(int i=k-1;i<m;i++){
                        (*sk)=(*sk)+pow(a[i][k-1], 2);
                }
                //CREA PIVOT DE CADA FILA
                pivot=a[k-1][k-1];
                //CREA SK
                if(pivot<=0) (*sk)=sqrt((*sk));
                else (*sk)=-sqrt((*sk));
                if ((*sk)<tol && (*sk)>-tol){
                    return -1;
                    break;
                }
                //CREA TK
                tau=((*sk)-pivot)/(*sk);
                //printf("\n\nTk=%g",tau);
                //CREA UK
                for (int i=0; i<=k-2; i++){
                    uk[i]=0;
                }
                uk[k-1]=1;
                //printf("\n\nU%d\n", k);
                for(int i=k;i<m;i++){
                    uk[i]=a[i][k-1]/(pivot-(*sk));
                    //printf("\n%g\t\n", uk[i]);
                }
                //Rellena el primer vector AK_i con el que genera Uk
                for(int i=k-1;i<m;i++){
                    a1[i]=a[i][k-1];
                }
                //Fa la primera columna amb els 0's canviats pel vector uk
                a[k-1][k-1]=a[k-1][k-1]-tau*prod_vectors(m, uk, a1)*uk[k-1];
                for(int j=k;j<m;j++){
                    a[j][k-1]=uk[j];
                }
                //RELLENAR la matriz AK
                cols=k;
                while(cols<n){
                    //GUARDA AK_i
                    for(int i=k-1;i<m;i++){
                        a1[i]=a[i][cols];
                    }
                    //ITERACION
                    for(int j=k-1;j<m;j++){
                        a[j][cols]=a[j][cols]-tau*prod_vectors(m, uk, a1)*uk[j];
                    }
                /*printf("\n\nmatriu A iteracio %d\n\n", cols+1);
                for(int i=0;i<m;i++){
                        printf("\n\n");
                    for(int j=0;j<n;j++){
                        printf("%.3f\t",a[i][j]);
                    }
                }*/
                cols++;
                }
            k++;
            }
            //RECUPERA MATRIZ A (sin ampliar)
            for(int i=0;i<m;i++){
                for(int j=0;j<n-1;j++){
                    matrix[i][j]=a[i][j];
                }
            }
            if ((*sk)<tol && (*sk)>-tol)return -1;
            /*printf("\n\n\n\nMatriu R diagonal superior\n");
            for(int i=0;i<m;i++){
                    printf("\n\n");
                for(int j=0;j<n;j++){
                    printf("%g\t",a[i][j]);
                }
            }

            printf("\n\n\n\nMatriu R diagonal superior(no ampliada)\n");
            for(int i=0;i<m;i++){
                    printf("\n\n");
                for(int j=0;j<n-1;j++){
                    printf("%g\t",matrix[i][j]);
                }
            }*/
            //RECUPERA VECTOR B
            //printf("\n\n\n\nVector Q^T*b\n\n");
                for(int j=0;j<m;j++){
                    b[j]=a[j][n-1];
                    //printf("%g\t",b[j]);
                }
            *tk=tau;
            //printf("\n\n\n\nTau=%g\n\n",*tk);
        free(a1);
        free(uk);
        esborra_matriu(a,m);
        return 0;
        }
    }
}

void matriu_dinamica(double ***m, int files, int cols)
{
    *m=malloc(files * sizeof(double *));
    for(int i=0;i<files;i++){
        (*m)[i]=malloc(cols * sizeof(double));
    }
}

void esborra_matriu(double **m, int files)
{
    for(int i=0;i<files;i++){
        free(m[i]);
    }
}

void vector_dinamic(double **vector, int size)
{
    *vector=malloc(size*sizeof(double));
}

double prod_vectors (int size, double *vf, double *vc)
{
    double r=0;
    for(int z=0;z<size;z++){
        r=r+((vf)[z]*((vc)[z]));
    }
    return r;
}

double **amplia_matriu(int f, int c, double **m, double *b)
{

    double **mAmpliada;
    matriu_dinamica(&mAmpliada,f, c+1);
    for(int j=0;j<f;j++){
        for(int z=0;z<c;z++){
        mAmpliada[j][z]=m[j][z];
        }
    }
    for(int z=0;z<f;z++){
        mAmpliada[z][c]=b[z];
    }
return mAmpliada;
}

double **prod_matrius (int f, int fc, int c, double **a, double **b)
{
    double **m;
    matriu_dinamica(&m, f,c);
    for(int i=0;i<f;i++){
        for(int j=0;j<c;j++){
            for(int z=0;z<fc;z++){
            m[i][j]=m[i][j]+((a)[i][z]*(b)[z][j]);
            }
        }
    }
return m;
}

void system_solve (int n, double **a, double *x)
{
    int dim=0;
    while(dim<n-1){
        x[n-2-dim]=a[n-2-dim][n-1];
        if (dim>0){
            for (int j=0; j<dim; j++){
                x[n-2-dim]=x[n-2-dim]-(a[n-2-dim][n-2-j]*x[n-2-j]);
            }
        }
    x[n-2-dim]=x[n-2-dim]/a[n-2-dim][n-2-dim];//PROBAMOS AQUI
    dim++;
    }
    //SOLUCIO DEL SISTEMA
    printf("\n\nSOLUCIO SISTEMA:\n\n");
    for(int i=0;i<n-1;i++){
        printf("%g\t",x[i]);
    }
}





int fdf(int m, int n, double *x, double *f, double **df)
{
    if (m!=3||n!=3)
    {
        printf("function1():: Error en les dimensions de la funció F.\n");
        return -1;
    }
    #define X1 x[0]
    #define X2 x[1]
    #define X3 x[2]
    #define F1 f[0]
    #define F2 f[1]
    #define F3 f[2]

    F1 = X1+X2*X3;
    F2 = X2-X1*X1;
    F3 = cos(X3);
    df[0][0] = 1;
    df[0][1] = X3;
    df[0][2] = X2;

    df[1][0] = -2*X1;
    df[1][1] = 1;
    df[1][2] = 0;

    df[2][0] = 0;
    df[2][1] = 0;
    df[2][2] = -sin(X3);

    /*printf("\n\n\n\nMatriu DF diagonal superior\n");
    for(int i=0;i<m;i++){
            printf("\n\n");
        for(int j=0;j<n;j++){
            printf("%g\t",df[i][j]);
        }
    }*/

    #undef X1
    #undef X2
    #undef X3
    #undef F1
    #undef F2
    #undef F3
return 0;
}

int newton(int m, int n, double *x, int (*fdf)(int m, int n, double *x, double *f, double **df), int maxit, double tol, double *f, double **df)
{

    int iter=0;
    double norma=0, tolQR=1e-10, tauk=0, sk=0;
    double *yn, **r;

    vector_dinamic(&yn, n);

    while(iter<maxit)
    {
        norma=0;

        fdf(m,n,x,f,df);//Dona valor a f i a df cada iterat

/*printf("Results from newton:");
printf("\n\nF(X):\n");
for (int i = 0; i < m; i++) {
    printf("%f\n", f[i]);
}

printf("\n\nDF(X):\n\n");
for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
        printf("%f\t ", df[i][j]);
    }
    printf("\n\n");
}*/

        //calcular el producte de df^-1*f con QRsolve;
        QRsolve(m, n+1, df, yn, tolQR, f, &tauk, &sk);


/*printf("Results from newton:");
printf("\n\nF(X):\n");
for (int i = 0; i < m; i++) {
    printf("%f\n", f[i]);
}

printf("\n\nDF(X):\n\n");
for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
        printf("%f\t ", df[i][j]);
    }
    printf("\n\n");
}*/
        r=amplia_matriu(m, n, df, f);

/*printf("\n\nMATRIU r:\n\n");
for (int i = 0; i < m; i++) {
    for (int j = 0; j < n+1; j++) {
        printf("%f\t ", r[i][j]);
    }
    printf("\n\n");
}*/

        if (sk>=tolQR || sk<=-tolQR)
        {
            system_solve(n+1,r,yn);
        }
        else return -1;

        for(int i=0; i<n; i++)
        {
            x[i]=x[i]-yn[i];
        }

        for (int i=0; i<m; i++)
        {
            norma=norma+pow(f[i], 2);
        }
        norma=sqrt(norma);
        if(norma<tol)
        {
            printf("\n\n\n\nHa convergit en %d iterats\n\n", iter+1);
            printf("\n\nSOLUCIO: X\n\n");
            for(int i=0;i<n;i++)
            {
                printf("%g\t",x[i]);
            }
            printf("\n\n\n\n");
        return 0;
        }
    iter++;
    }

    if(iter==maxit)
    {
        printf("\n\nS'ha esgotat el nombre maxim d'iterats\n\n");
        return -1;
    }

    esborra_matriu(r, m);
    free(yn);
return -1;
}






void camp1(int n, double t, double *x, double *f, Prm *prm)
{

    double a=prm->a;

    f[0]=a*x[0]-x[1];
    f[1]=x[0]+a*x[1];
}

void solucioExacta(double t, double x0, double y0, double *xExacta, Prm *p)
{

    xExacta[0]=exp((p->a)*t)*(x0*cos(t)-y0*sin(t));
    xExacta[1]=exp((p->a)*t)*(x0*sin(t)+y0*cos(t));
}



int rk45(double *t, double *x,  double *h, double hmin, double hmax, double tol, int n, void (*camp)(int n, double t, double *x, double *f, Prm *prm), Prm *prm)
{

    //Coeficients RK45 (estatics)
    static double c[6] = {0., 1./4., 3./8., 12./13., 1., 1./2.};
    static double A[15] = {1./4., 3./32., 9./32., 1932./2197.,
        -7200./2197.,7296./2197., 439./216., -8., 3680./513.,
        -845./4104., -8./27., 2., -3544./2565., 1859./4104., -11./40.};
    static double bp[6] = {25./216., 0., 1408./2565., 2197./4104., -1./5.,0.};
    static double b[6] = {16./135., 0., 6656./12825., 28561./56430., -9./50., 2./55.};

    double ek;
    int s=6, aux;
    double x4[n], x5[n], ki[n], F4[n], F5[n], xki[n], Mki[n][s];

    while(1)
    {
        aux=0;
        ek=0;

        //Omplim el ki inicial que es f(tk,xk)
        camp(n,*t+c[0]*(*h),x,ki,prm);

        for(int j=0;j<n;j++)
        {
            xki[j]=x[j];
            Mki[j][0]=ki[j];
            F4[j]=0;
            F5[j]=0;
        }

        //Calculem els Ki's
        for(int i=1;i<s;i++)
        {
            if(i>0) aux+=i-1;
            for(int j=0;j<n;j++)
            {
                for(int l=0;l<i;l++)
                {
                    xki[j]+=(*h)*A[aux+l]*Mki[j][l];//Se actualiza la xk (xki)
                }
            }
            camp(n,*t+c[i]*(*h),xki,ki,prm);//Se actualiza ki para siguiente paso
            for(int j=0;j<n;j++)
            {
                Mki[j][i]=ki[j];
                xki[j]=x[j];
            }
        }

        //Calculem les F's
        for(int j=0;j<n;j++)
        {
            for(int i=0;i<s;i++)
            {
                F4[j]+=bp[i]*Mki[j][i];
                F5[j]+=b[i]*Mki[j][i];
            }
        }

        //Calculem les solucions d'ordre 4 i 5
        for(int j=0;j<n;j++)
        {
            x4[j]=x[j]+(*h)*F4[j];
            x5[j]=x[j]+(*h)*F5[j];
            //Calculem error
            ek+=(x5[j]-x4[j])*(x5[j]-x4[j]);
        }
        ek=sqrt(ek);

        //Verificar si l'error es acceptable
        if(ek<=tol)

        {
            //Acceptem x_p+1 el valor de x(tk+1)
            for(int i=0;i<n;i++)
            {
                x[i]=x5[i];
            }
            *t+=*h;

            //Nou pas h_k+1
            if(ek<tol/256) ek=tol/256;
            (*h)*=0.9*pow(tol/ek, 0.125);
            if(fabs(*h)>hmax)
            {
                *h=copysign(hmax,*h);
            }

            if(fabs(*h)<hmin)
            {
                printf("\nNo es pot calcular amb la precisió demanada, ja que hmin es massa gran\n");
                return -1;
            }
        break;
        }
        else
        {
            if(*h==hmin)
            {
                printf("\nNo es pot calcular amb la precisió demanada, ja que hmin es massa gran\n");
                return -1;
            }
            else if(fabs(*h/2)<hmin) *h=copysign(hmin,*h);
            else if(fabs(*h/2)>=hmin) *h=*h/2;
        }
    }
return 0;
}

int rk78(double *t, double *x, double *h, double hmin, double hmax, double tol, int n, void (*camp)(int n, double t, double *x, double *f, Prm *prm), Prm *prm)
{

    // Coeficients RK78 (estatics)
    static double c[13] = {0., 2./27., 1./9., 1./6., 5./12., 1./2., 5./6., 1./6., 2./3., 1./3., 1., 0., 1.};
    static double A[78] = {
        2./27.,
        1./36., 1./12.,
        1./24., 0., 1./8.,
        5./12., 0., -25./16., 25./16.,
        1./20., 0., 0., 1./4., 1./5.,
        -25./108., 0., 0., 125./108., -65./27., 125./54.,
        31./300., 0., 0., 0., 61./225., -2./9., 13./900.,
        2., 0., 0., -53./6., 704./45., -107./9., 67./90., 3.,
        -91./108., 0., 0., 23./108., -976./135., 311./54., -19./60., 17./6., -1./12.,
        2383./4100., 0., 0., -341./164., 4496./1025., -301./82., 2133./4100., 45./82., 45./164., 18./41.,
        3./205., 0., 0., 0., 0., -6./41., -3./205., -3./41.,3./41., 6./41., 0.,
        -1777./4100., 0., 0., -341./164., 4496./1025., -289./82., 2193./4100., 51./82., 33./164., 12./41., 0., 1.};

    static double bp[13] = {41./840., 0., 0., 0., 0., 34./105., 9./35., 9./35., 9./280., 9./280., 41./840., 0., 0.};
    static double b[13] = {0., 0., 0., 0., 0., 34./105., 9./35., 9./35., 9./280., 9./280., 0., 41./840., 41./840.};

    double ek;
    int s=13, aux;
    double x7[n], x8[n], ki[n], F7[n], F8[n], xki[n], Mki[n][s];

    while(1)
    {
        aux=0;
        ek=0;

        // Omplim el ki inicial que es f(tk,xk)
        camp(n,*t+c[0]*(*h),x,ki,prm);

        for(int j=0;j<n;j++)
        {
            xki[j]=x[j];
            Mki[j][0]=ki[j];
            F7[j]=0;
            F8[j]=0;
        }

        // Calculem els Ki's
        for(int i=1;i<s;i++)
        {
            if(i>0) aux+=i-1;
            for(int j=0;j<n;j++)
            {
                for(int l=0;l<i;l++)
                {
                    xki[j]+=(*h)*A[aux+l]*Mki[j][l]; // Se actualiza la xk (xki)
                }
            }
            camp(n,*t+c[i]*(*h),xki,ki,prm); // Se actualiza ki para siguiente paso
            for(int j=0;j<n;j++)
            {
                Mki[j][i]=ki[j];
                xki[j]=x[j];
            }
        }

        // Calculem les F's
        for(int j=0;j<n;j++)
        {
            for(int i=0;i<s;i++)
            {
                F7[j]+=bp[i]*Mki[j][i];
                F8[j]+=b[i]*Mki[j][i];
            }
        }

        // Calculem les solucions d'ordre 7 i 8
        for(int j=0;j<n;j++)
        {
            x7[j]=x[j]+(*h)*F7[j];
            x8[j]=x[j]+(*h)*F8[j];
            // Calculem error
            ek+=(x8[j]-x7[j])*(x8[j]-x7[j]);
        }
        ek=sqrt(ek);

        // Verificar si l'error es acceptable
        if(ek<=tol)
        {
            // Acceptem x_p+1 el valor de x(tk+1)
            for(int i=0;i<n;i++)
            {
                x[i]=x8[i];
            }
            *t+=*h;

            // Nou pas h_k+1
            if(ek<tol/256) ek=tol/256;
            (*h)*=0.9*pow(tol/ek, 0.125);
            if(fabs(*h)>hmax)
            {
                *h=copysign(hmax,*h);
            }

            if(fabs(*h)<hmin)
            {
                printf("\nNo es pot calcular amb la precisió demanada, ja que hmin es massa gran\n");
                return -1;
            }
        break;
        }
        else
        {
            if(*h==hmin)
            {
                printf("\nNo es pot calcular amb la precisió demanada, ja que hmin es massa gran\n");
                return -1;
            }
            else if(fabs(*h/2)<hmin) *h=copysign(hmin,*h);
            else if(fabs(*h/2)>=hmin) *h=*h/2;
        }
    }
return 0;
}


void lorenz_system(int n, double t, double *x, double *f, Prm *prm)
{
    double sigma = prm->sigma;
    double rho = prm->rho;
    double beta = prm->beta;

    f[0] = sigma * (x[1] - x[0]);
    f[1] = x[0] * (rho - x[2]) - x[1];
    f[2] = x[0] * x[1] - beta * x[2];
}

void variacionales4(int n, double t, double *x, double *f, Prm *prm)
{
    double alpha = prm->a;
    double r2 = x[0] * x[0] + x[1] * x[1];

    if(n==6)
    {
        for (int k = 0; k < 6; k++)
        {
            f[k] = 0;
        }

        // EDO principales
        f[0] = alpha * (1 - r2) * x[0] - x[1];
        f[1] = x[0] + alpha * (1 - r2) * x[1];

        // Jacobiano del sistema
        double J[2][2];
        J[0][0] = alpha * (1 - r2) - 2 * alpha * x[0] * x[0];
        J[0][1] = -1 - 2 * alpha * x[0] * x[1];
        J[1][0] = 1 - 2 * alpha * x[0] * x[1];
        J[1][1] = alpha * (1 - r2) - 2 * alpha * x[1] * x[1];

        // Ecuaciones variacionales (A)
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    f[2 + (j) * 2 + (i)] += J[i][k]*x[2 + (j) * 2 + (k)];
                }
            }
        }
    }

}

int flux(double *t, double *x, double *h, double T, double hmin, double hmax, double tol, int npasmax, int n, void (*camp)(int n, double t, double *x, double *f, Prm *prm), Prm *prm, int identif)
{
    double t_end = *t + T;
    int passos = 0;
    double *f;

    f=malloc(n*sizeof(double));


    if((T > 0 && *h < 0) || (T < 0 && *h > 0)){*h = -*h;}


    FILE *output1=fopen("lorenz_data1.dat", "w");
    if (output1==NULL) {
        perror("No es pot crear l'arxiu");
        return -1;
    }

    FILE *output2 = fopen("pr4ex3_results.dat", "w");
    if (output2==NULL){
        printf("Error: No se pudo abrir el archivo para escribir resultados\n");
        return -1;
    }
    //fprintf(output2, "t x y A11 A12 A21 A22\n");


    // Guardar la órbita en un archivo para graficar
    FILE *output3 = fopen("orbHalo.txt", "w");
    if (!output3) {
        fprintf(stderr, "Error: No se pudo crear el archivo orbHalo.txt\n");
        return -1;
    }

    while((T > 0 && *t < t_end) || (T < 0 && *t > t_end))
    {
        if(passos >= npasmax)
        {
            printf("\nS'ha esgotat el nombre maxim de passosFlux\n");
            return -10;
        }

        if(fabs(*t + *h) > fabs(t_end))
        {
            *h = t_end - *t;
            *t=t_end;
        }

        if(identif==42){fprintf(output1, "%f %f %f\n", x[0], x[1], x[2]);}
        else if(identif==43){fprintf(output2, "%f %f %f %f %f %f %f\n", *t, x[0], x[1], x[2], x[3], x[4], x[5]);}
        else if(identif==55){fprintf(output3, "%f %.16f %.5f %.18f %.5f %.16f %.5f\n", *t, x[0], x[1], x[2], x[3], x[4], x[5]);}


        //Fer anar rk78
        if(rk78(t, x, h, hmin, hmax, tol, n, camp, prm)!=0)
        {
            fprintf(stderr,"Error: No se pudo alcanzar la precision con hmin=%g en t%d=%g\n",hmin,passos,*t);
            break;
        }

    passos++;
    }

    if(identif==42){fclose(output1);}
    else if(identif==43){fclose(output2);}
    else if(identif==55){fclose(output3);}
    free(f);

    return 0;
}



int Ham(double *x, double *H, Prm *prm)
{
    *H = -cos(x[0]) + 0.5 * x[1] * x[1];
    return 0;
}

void pendolvariacional(int n, double t, double *x, double *f, Prm *prm)
{
   f[0]=x[1];
   f[1]=-sin(x[0]);
   if (n==6) {
      f[2]=x[3];
      f[3]=-x[2]*cos(x[0]);
      f[4]=x[5];
      f[5]=-x[4]*cos(x[0]);
   }
}



int opham_fdf(int m, double hh, double *xx, double *cg, double *ff, double **df,double h0rk, double hminrk, double hmaxrk, double tolrk, int npasmxrk,int (*ham)(double *x, double *H, Prm *prm),void (*camp)(int n, double t, double *x, double *f, Prm *prm),Prm *prm)
{

    int n = 2 * m;
    double trk = 0.0;
    double t=0.0; //Es la t de camp
    double *x = &xx[1];
    double *DH = malloc(n*sizeof(double));
    double T = xx[0];

    double H=0;
    double gx = cg[n];
    double *phi_T = malloc(n * sizeof(double));
    double *phi_TA = malloc((n+n*n) * sizeof(double));
    double *Dgx = malloc(n*sizeof(double));
    double *fphi_T = malloc(n*sizeof(double));
    double **Dphi_T;
    double *DHAux = malloc(m*sizeof(double));

    matriu_dinamica(&Dphi_T, n, n);

    /*printf("\n\n");
    printf("XX\n");
    for (int i = 0; i < n+1; i++) {
        printf("xx[%d]=%f\t",i,xx[i]);
    }
    printf("\n\n");

    printf("\n\n");
    printf("X\n");
    for (int i = 0; i < n; i++) {
        printf("x[%d]=%f\t",i,x[i]);
    }
    printf("\n\n");*/

    //calculamos H(x)=H
    ham(x, &H, prm);
    /*printf("\n\n");
    printf("H=%f\t",H);
    printf("\n\n");*/

    // Calculamos g(x)
    for (int i = 0; i < n; i++) {
        gx += cg[i]*x[i];
    }
    /*printf("\n\n");
    printf("gx=%f\t",gx);
    printf("\n\n");*/

    // Calculamos phi_T(x) y Dphi_T(x)
    for (int i = 0; i < n; i++) {
        phi_TA[i] = x[i];
    }
    /*printf("\n\n");
    printf("X(0)\n");
    for (int i = 0; i < n; i++) {
        printf("phi_T0[%d]=%f\t",i,phi_TA[i]);
    }
    printf("\n\n");*/

    //GUARDAR LA ID
    int auxphi=0;
    for (int i = 0; i < n*n; i++)
    {
        if(i==auxphi*n+auxphi)
        {
            phi_TA[i+n]=1;
            auxphi++;
        }
        else phi_TA[i+n]=0;
    }
    /*printf("\n\nDphi_T0\n\n");
    for (int i = 0; i < n*n; i++)
    {
        if(i==n || i==2*n || i==3*n || i==4*n || i==5*n) printf("\n");
        printf("%f\n",phi_TA[i+n]);
    }
    printf("\n\n");*/

    // Calculamos phi_T(x) y Dphi_T(x)
    flux(&trk, phi_TA, &h0rk, T, hminrk, hmaxrk, tolrk, npasmxrk, n+n*n, camp, prm, 50);

    //phi_TA -> phi_T
    for (int i = 0; i < n; i++) {
        phi_T[i] = phi_TA[i];
    }
    /*printf("\n\n");
    printf("FLUX\n");
    for (int i = 0; i < n; i++) {
        printf("phi_T[%d]=%f\t",i,phi_TA[i]);
    }
    printf("\n\n");

    printf("\n\nDphi_T\n\n");
    for (int i = 0; i < n*n; i++)
    {
        if(i==n || i==2*n || i==3*n || i==4*n || i==5*n) printf("\n");
        printf("%.2f\n",phi_TA[i+n]);
    }
    printf("\n\n");*/

    // Calculamos F(X)
    ff[0] = H - hh;
    ff[1] = gx;
    for (int i = 0; i < n; i++) {
        ff[i + 2] = phi_T[i] - x[i];
    }
    /*printf("\n\n");
    printf("F(x)\n");
    for (int i = 0; i < n; i++) {
        printf("ff[%d]=%f\t",i,ff[i]);
    }
    printf("\n\n");*/

    // Calculamos DF(X) if df != NULL
    if (df != NULL)
    {

        // Calculamos DH(x)
        camp(n, t, x, DH, prm);

        /*printf("\n\nDH/f(x)\n");
        for (int i = 0; i < n; i++)
        {
            printf("DH[%d]=%f\t",i,DH[i]);
        }
        printf("\n\n");*/

        for(int i=0; i<m; i++)
        {
            DHAux[i]=DH[i];
            DH[i]=-DH[i+m];
            DH[i+m]=DHAux[i];
        }

        /*printf("\n\nDH BUENA\n");
        for (int i = 0; i < n; i++)
        {
            printf("DH[%d]=%f\t",i,DH[i]);
        }
        printf("\n\n");*/

        // Calculamos Dg(x) esta es constante y no cambia ya que no depende de x
        for (int i = 0; i < n; i++)
        {
            Dgx[i] = cg[i];
        }

        /*printf("\n\nDgx\n");
        for (int i = 0; i < n; i++)
        {
            printf("Dgx[%d]=%f\t",i,Dgx[i]);
        }
        printf("\n\n");*/

        // Calculamos f(phi_T(x))
        camp(n,t,phi_T,fphi_T,prm);

        /*printf("\n\nf(phi_T)\n");
        for (int i = 0; i < n; i++)
        {
            printf("fphi_T[%d]=%f\t",i,fphi_T[i]);
        }
        printf("\n\n");*/


        // Calculamos Dphi_T(x)-I
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if(i==j)
                {
                    Dphi_T[j][i] = phi_TA[n+i*n+j]-1;
                }
                else Dphi_T[j][i] = phi_TA[n+i*n+j];
            }
        }

        /*printf("\n\nMATRIU Dphi_T-Id\n");
        for (int i = 0; i < n; i++)
        {
            printf("\n\n");
            for (int j = 0; j < n; j++)
            {
                printf("%.4f\t\t",Dphi_T[i][j]);
            }
        }
        printf("\n\n");*/



        // Calculamos df (DF(x))

        //METEMOS DH
        df[0][0]=0;
        for (int j = 0; j < n; j++)
        {
            df[0][j+1] = DH[j];
        }

        //METEMOS Dgx
        df[1][0]=0;
        for (int j = 0; j < n; j++)
        {
            df[1][j+1] = Dgx[j];
        }

        //METEMOS fphi_T(x) y Dphi_T(x)-I
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n+1; j++)
            {
                if(j==0)
                {
                    df[i+2][j] = fphi_T[i];
                }
                else df[i+2][j] = Dphi_T[i][j-1];
            }
        }
    }

        /*printf("\n\nMATRIU DF\n");
        for (int i = 0; i < n+2; i++)
        {
            printf("\n\n");
            for (int j = 0; j < n+1; j++)
            {
                printf("%.4f\t\t",df[i][j]);
            }
        }
        printf("\n\n");*/

    free(DH);
    free(phi_T);
    free(phi_TA);
    free(Dgx);
    free(fphi_T);
    free(DHAux);
    esborra_matriu(Dphi_T,n);

return 0;
}



int newton1(int (*opham_fdf)(int m, double hh, double *xx, double *cg, double *ff, double **df, double h0rk, double hminrk, double hmaxrk, double tolrk, int npasmxrk, int (*ham)(double *x, double *H, Prm *prm), void (*camp)(int n, double t, double *x, double *f, Prm *prm), Prm *prm),int m, double hh, double *xx, double *cg, double *ff, double **df, double h0rk, double hminrk, double hmaxrk, double tolrk, int npasmxrk,int (*ham)(double *x, double *H, Prm *prm),void (*camp)(int n, double t, double *x, double *f, Prm *prm),int maxit, double tol, Prm *prm)
{

    int n=2*m;

    int iter = 0;
    double norma = 0, tolQR = 1e-20, tauk = 0, sk = 0;
    double norm_dx=0;
    double *yn, **r;

    // Asignación dinámica del vector auxiliar
    vector_dinamic(&yn, n+1);

    while (iter < maxit)
    {
        norm_dx=0;
        norma = 0;

        // Evaluación de la función y su jacobiano
        opham_fdf(m, hh, xx, cg, ff, df, h0rk, hminrk, hmaxrk, tolrk, npasmxrk, ham, camp, prm);

        // Resolver el sistema lineal usando QRsolve
        QRsolve(n+2, n+2, df, yn, tolQR, ff, &tauk, &sk);

        // Ampliar la matriz jacobiana con los valores de la función
        r = amplia_matriu(n+2, n+1, df, ff);

        // Verificar la singularidad
        if (sk >= tolQR || sk <= -tolQR)
        {
            system_solve(n+2, r, yn);
        }
        else return -10;

        // Actualizar las soluciones
        for (int i = 0; i < n+1; i++) {
            xx[i] -= yn[i];
        }

        // Calcular la norma de F(X)
        for (int i = 0; i < n+2; i++) {
            norma += ff[i] * ff[i];
        }
        norma = sqrt(norma);

        // Calcular la norma de X
        for (int i = 0; i < n + 1; i++) {
            norm_dx += yn[i] * yn[i];
        }
        norm_dx = sqrt(norm_dx);

        // Imprimir iteración, norma de F(X) y norma de X
        printf("newton: iter %d nf %.8f nc %.8f", iter, norma, norm_dx);

        // Verificar criterio de convergencia
        if (norma < tol) {return iter;}

    iter++;
    }

    // Si no converge, informar y liberar memoria
    if (iter == maxit) {
        printf("\n\nSe alcanzo el numero maximo de iteraciones sin convergencia.\n\n");
    }

    esborra_matriu(r, n+2);
    free(yn);
    return -1;
}



int opham(int m, double hh, double *xx, double tol, int maxit, double *cg,double h0rk, double hminrk, double hmaxrk, double tolrk, int npasmxrk,int (*ham)(double *x, double *H, Prm *prm),void (*camp)(int n, double t, double *x, double *f, Prm *prm),Prm *prm)
{

    int n = 2 * m;
    double *ff , **df;

    vector_dinamic(&ff,n+2);
    matriu_dinamica(&df, n+2, n+1);

    int iter = newton1(opham_fdf, m, hh, xx, cg, ff, df, h0rk, hminrk, hmaxrk, tolrk, npasmxrk, ham, camp, maxit, tol, prm);

    if (iter >= 0) {
        return iter;
    } else {
        return -1;
    }

    esborra_matriu(df,n+2);
    free(ff);
}


int rtbp_h (double *x, double *H, Prm *prm)
{
   double mu=prm->mu, xmmu=X-mu, xmmup1=xmmu+1,
	  r12=SQR(xmmu)+SQR(Y)+SQR(Z),
	  r22=SQR(xmmup1)+SQR(Y)+SQR(Z),
	  r1=sqrt(r12), r2=sqrt(r22), /*r13=r1*r12, r23=r2*r22,*/
	  p1=(1-mu)/r1, p2=mu/r2;
   /*if (dmu!=NULL) *dmu=1/r1-1/r2-xmmu*p1/r12-xmmup1*p2/r22;*/
   *H=.5*(SQR(PX)+SQR(PY)+SQR(PZ))+Y*PX-X*PY-p1-p2;
   return 0;
}


void rtbp(int n, double t, double *x, double *f, Prm *prm)
{
   double r1x, r2x, r1x2, r2x2, r1ypz, r12, r22, r13, r23, p13, p23,
	  p123, p15, p25, p125, p15x, p25x, p125x,
	  dxpx, dypx, dzpx, dypy,
	  dzpy, dzpz;
    double MU=prm->mu;
   int j;
/* Equacions */
   f[0]=PX+Y; f[1]=PY-X; f[2]=PZ;
   r1x=X-MU; r2x=r1x+1; r1x2=r1x*r1x; r2x2=r2x*r2x;
   r1ypz=Y*Y+Z*Z;
   r12=r1x2+r1ypz; r22=r2x2+r1ypz;

   r13=r12*sqrt(r12); r23=r22*sqrt(r22);
   p13=(1-MU)/r13; p23=MU/r23;

   f[3]=PY-(p13*r1x+p23*r2x);
   p123=p13+p23;
   f[4]=-PX-Y*p123;
   f[5]=-Z*p123;
   if (n>=NV1) {
   /* Variacionals primeres */
      p15=p13/r12; p25=p23/r22;
      p125=p15+p25;
      p15x=p15*r1x; p25x=p25*r2x;
      p125x=p15x+p25x;
      dxpx=-p123+3*(p15x*r1x+p25x*r2x);
      dzpx=dypx=3*p125x;
      dypx*=Y; dzpx*=Z;
      dzpy=dypy=3*Y*p125;
      dypy=Y*dypy-p123;
      dzpy*=Z;
      dzpz=-p123+3*Z*Z*p125;
  /* Derivada resp var dep del camp per der flux resp c.i. (part de var dep) */
      for (j=0; j<N; j++) {
#define VARF(i,j) f[(1+(j))*N+(i)]
#define VAR(i,j) x[(1+(j))*N+(i)]
	 VARF(0,j)= VAR(1,j)+VAR(3,j);
	 VARF(1,j)=-VAR(0,j)+VAR(4,j);
	 VARF(2,j)= VAR(5,j);
	 VARF(3,j)=dxpx*VAR(0,j)+dypx*VAR(1,j)+dzpx*VAR(2,j)
	           +VAR(4,j);
	 VARF(4,j)=dypx*VAR(0,j)+dypy*VAR(1,j)+dzpy*VAR(2,j)
	           -VAR(3,j);
	 VARF(5,j)=dzpx*VAR(0,j)+dzpy*VAR(1,j)+dzpz*VAR(2,j);
#undef VAR
#undef VARF
      }
   }
}
