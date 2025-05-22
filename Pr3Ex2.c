#include "FuncionsQR.h"


int Pr3Ex2(void){

    int n=2, passos=0;
    double hmin=0.001, hmax=0.01, tol=1e-8, t, h=0.01, x0=1.0, y0=0.0;
    double *x, *f, *solExc;

    Prm *P1=malloc(sizeof(Prm));
    Prm *P2=malloc(sizeof(Prm));
    Prm *P3=malloc(sizeof(Prm));
    x=malloc(n*sizeof(double));
    f=malloc(n*sizeof(double));
    solExc=malloc(n*sizeof(double));

    //Guardem valor dels parametres
    P1->a=0; P1->tMax=2*M_PI;
    P2->a=0.1; P2->tMax=2*M_PI;
    P3->a=0.1; P3->tMax=-20*M_PI;

    //Establim condicions inicials
    x[0]=x0;
    x[1]=y0;

    t=0;

    //Calcul de f en x0,y0
    camp1(n,t,x,f,P1);

    //Arxiu per guardar els resultats
    FILE *output=fopen("Resultat_rk78.txt","w");
    if(output==NULL)
    {
        perror("No es pot crear l'arxiu");
        return 1;
    }

    while (t<P1->tMax)
    {

        if(passos==50000) printf("\nS'ha esgotat el nombre maxim de passos\n");
        //Calcular la solució exacta
        solucioExacta(t,x0,y0,solExc,P1);

        //Calcular la solució exacta
        //solucioExacta(t,x0,y0,solExc,P2);

        //Calcular la solució exacta
        //solucioExacta(t,x0,y0,solExc,P3);

        //Escriure els resultats de cada pas en l'arxiu
        fprintf(output,"t%d = %g, x%d = %g, y%d = %g, x = %g, y = %g, error_x = %e, error_y = %e\n",passos,t,passos,x[0],passos,x[1],solExc[0],solExc[1],fabs(x[0]-solExc[0]),fabs(x[1]-solExc[1]));

        //Fer anar rk45
        if(rk78(&t,x,&h,hmin,hmax,tol,n,camp1,P1)!=0)
        {
            fprintf(stderr,"Error: No se pudo alcanzar la precision con hmin=%g en t%d=%g\n",hmin,passos,t);
            break;
        }

        //Fer anar rk45
        /*if(rk45(&t,x,&h,hmin,hmax,tol,n,camp,P2)!=0)
        {
            fprintf(stderr,"Error: No se pudo alcanzar la precision con hmin=%g en t%d=%g\n",hmin,k,t);
            break;
        }*/


        //Fer anar rk45
        /*if(rk45(&t,x,&h,hmin,hmax,tol,n,camp,P3)!=0)
        {
            fprintf(stderr,"Error: No se pudo alcanzar la precision con hmin=%g en t%d=%g\n",hmin,k,t);
            break;
        }*/

    passos++;
    }

    fclose(output);

    printf("Resultados guardados en Resultat_rk78.txt\n");

return 0;
}

