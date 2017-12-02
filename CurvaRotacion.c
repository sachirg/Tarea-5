#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

int number_lines(void);

/*-----------------------------------------------------
FUNCIONES
---------------------------------------------------------*/

void Vcircular(double *_vcircular ,double *Radio, double Mb, double Md, double Mh);
//Likelihood de chi cuadrado
double likelihood(double *y_obs, double *y_model);
//load_data
void load_data(double *radios, double *velocidadesObservadas);
//Punto aleatorio dist. Uniforme
double drand();



//obtener el maximo de un arreglo
int max(double *vector, int size);
//Imprimir condiciones iniciales
void imprimirIniciales(double *_vector1, double *_vector2);
//Imprimir arreglo en vector.txt
void imprimirVector(double *_vector);
//Imprimir arreglo en vector2.txt
void imprimirVector2(double *_vector);



/*-----------------------------------------------------
CONSTANTES
---------------------------------------------------------*/

double b_b=0.2497;  // [KPC]
double b_d=5.16;    // [KPC] 
double a_d=0.3105;  // [KPC]
double a_h=64.3;    // [KPC]
double mGal;        // [KPC]
double parsec;      // [KPC]

/*-----------------------------------------------------
PUNTEROS
---------------------------------------------------------*/

char *name;
double *radios;
double *velocidadesObservadas;
double *M_d_walk;
double *M_b_walk;
double *M_h_walk;
double *v_0,*v_1,*l_walk;

/*-----------------------------------------------------
ESTIMADORES
---------------------------------------------------------*/
double l_0,l_1, alpha,beta,Mb_prime,Md_prime,Mh_prime;
int Niteraciones;
int imax;
int nDatos;


int main(void)
{   


    parsec=1/pow(3.087,16);
    mGal==pow(2.325,7);

    nDatos=number_lines()-1;

    Niteraciones=50000;

    double n=1;

    M_b_walk=malloc(Niteraciones*sizeof(double));
    M_d_walk=malloc(Niteraciones*sizeof(double));
    M_h_walk=malloc(Niteraciones*sizeof(double));


    M_b_walk[0]=drand();
    M_d_walk[0]=drand();
    M_h_walk[0]=drand();

    radios=malloc(nDatos*sizeof(double));
    velocidadesObservadas=malloc(nDatos*sizeof(double));
    v_0=malloc(nDatos*sizeof(double));
    v_1=malloc(nDatos*sizeof(double));


    l_walk=malloc(Niteraciones*sizeof(double));

    load_data(radios,velocidadesObservadas);

    Vcircular(v_0, radios, M_b_walk[0], M_d_walk[0], M_h_walk[0]);

 
 
    int j;
    for (j=1;j<Niteraciones;j++)
    {
        Mb_prime=M_b_walk[j-1]+drand();
        Md_prime=M_d_walk[j-1]+drand();
        Mh_prime=M_h_walk[j-1]+drand();        

        free(v_0);
        free(v_1);
        v_0=malloc(nDatos*sizeof(double));
        v_1=malloc(nDatos*sizeof(double));


        Vcircular(v_0,radios, M_b_walk[j-1], M_d_walk[j-1], M_h_walk[j-1]);
        Vcircular(v_1,radios, Mb_prime, Md_prime, Mh_prime);

        l_0=likelihood(velocidadesObservadas,v_0);
        l_walk[j]=l_0;
        l_1=likelihood(velocidadesObservadas,v_1);

        
        alpha=l_1/l_0;

        if(alpha>=1.0)
        {

            M_b_walk[j+1]=Mb_prime+(2*drand()-1);
            M_d_walk[j+1]=Md_prime+(2*drand()-1);
            M_h_walk[j+1]=Mh_prime+(2*drand()-1);        
        }
        else
        {
            beta=drand();
            if(beta<alpha)
            {
                M_b_walk[j+1]=Mb_prime;
                M_d_walk[j+1]=Md_prime;
                M_h_walk[j+1]=Mh_prime; 
            }
            else
            {
                M_b_walk[j+1]=M_b_walk[j-1]+(2*drand()-1);
                M_d_walk[j+1]=M_d_walk[j-1]+(2*drand()-1);
                M_h_walk[j+1]=M_h_walk[j-1]+(2*drand()-1); 
            }
        }

    }

    int mejori=max(l_walk,Niteraciones);
    Vcircular(v_0,radios, M_b_walk[mejori], M_d_walk[mejori], M_h_walk[mejori]);
    imprimirVector(v_0);

    printf("Los mejores parametros son: \n Mb:%e\n Md:%e\n Mh:%e \n",M_b_walk[mejori], M_d_walk[mejori], M_h_walk[mejori]);

    return 0;

}


void load_data(double *x, double *y)
{
    char *delimiter=" ";

    FILE* f = fopen("RadialVelocities.dat", "r");

    char palabra[100];
    double number;
    int i=0;

    fscanf(f, "%s",palabra);
    fscanf(f, "%s",palabra);

    while(fscanf(f, "%lf",&number)==1)
    {
        if(i%2==0)
        {
            radios[i/2] = number;
        }
        else
        {

            velocidadesObservadas[i/2] = number;
        }
        i++;

    }


    fclose(f);

}


void imprimirIniciales(double *_vector1, double *_vector2)
{
    FILE *fout;
    fout=fopen("iniciales.txt","w+");
    for(int i=0;i<nDatos-1;i++)
    {
            fprintf(fout,"%lf %lf\n",_vector1[i],_vector2[i]);
    }
    
}

void imprimirVector(double *_vector)
{
    FILE *fout;
    fout=fopen("daticos.txt","w+");
    for(int i=0;i<nDatos-1;i++)
    {
            fprintf(fout,"%lf\n",_vector[i]);
    }
    
}

void Vcircular(double *_vcircular ,double *_Radio, double _Mb, double _Md, double _Mh)
{

    int i;


    for(i=0;i<nDatos-1;i++)
    {
        double t1,t2,t3,t1arriba,t1abajo,t2arriba,t2abajo,t3arriba,t3abajo;
        
        t1arriba=sqrt(_Mb)*_Radio[i];
        t1abajo=pow(_Radio[i],2) + pow(b_b,2);
        t1abajo=pow(t1abajo,3);
        t1abajo=sqrt(t1abajo);
        t1abajo=sqrt(t1abajo);

        t2arriba=sqrt(_Md)*_Radio[i];
        t2abajo=pow(_Radio[i],2) + pow((b_d+a_d),2);
        t2abajo=pow(t2abajo,3);
        t2abajo=sqrt(t2abajo);
        t2abajo=sqrt(t2abajo);

        t3arriba=sqrt(_Mh);
        t3abajo=pow(_Radio[i],2) + pow(a_h,2);
        t3abajo=sqrt(t3abajo);
        t3abajo=sqrt(t3abajo);

        t1=t1arriba/t1abajo;
        t2=t2arriba/t2abajo;
        t3=t3arriba/t3abajo;

        _vcircular[i]=t2+t1+t3;
    }
}


double likelihood(double *y_obs, double *y_model)
{

    int i =0;
    double chi_squared=0;
    for (i=0;i<nDatos-1;i++)
    {
        chi_squared=chi_squared+pow((y_obs[i]-y_model[i]),2);
    }
    chi_squared=(1.0/2.0)*chi_squared;
    return exp(-(chi_squared)/nDatos);
}


double drand()
{
  return (rand()+1.0)/(RAND_MAX+1.0);
}


//Contar numero de lineas
int number_lines(void)
{
    int lines = 0, ch;

    FILE *file = fopen("RadialVelocities.dat", "r");
    
    while(!feof(file))
    {

        ch = fgetc(file);
        if(ch == '\n')
        {
            lines++;

        }
    }

    return lines;
}



int max(double *vector, int size)
{
        double max = vector[0];
        int max_i=0;
        int i;
        for(i = 1; i < size; ++i)
        {
            if(vector[i] > max)
            {
                max = vector[i];
                max_i = i;
            }
        }      
        return max_i;
 }
