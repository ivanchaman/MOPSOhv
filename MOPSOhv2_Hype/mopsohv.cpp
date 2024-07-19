/*Acotado*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "randomlib.h"
#include "indicadorHype.h"

#define popsize      100   /* numero de particulas en la poblacion */
#define archive_size 100   /* tamaño del archivo */
#define verbose      1
#define printevery   10
#define PI 3.14159265358979
#define RMIN  -1.18e-38
#define acotado 1

int function;       /* codigo de la funcion objetivo */
int maxgen;         /* numero maximo de generaciones */
int optimization;   /* tipo de optimizacion, 0 minimizacion, 1 maximizacion*/
int maxfun;         /* numero maximo de funciones objetivo */
int maxvar;         /* numero maximo de variables */

double bound = 11.0;

//double  PI;
double **archiveVar; /* valores de las variables de las particulas en el archivo*/
double **archiveFit; /* valores de aptitud de las particulas en el archivo */

double **popVar;	 /* valores de las variables de las particulas en la poblacion*/
double **popFit;     /* valores de aptitud de las particulas en la poblacion*/
double **pbestsVar;	 /* valor pBest de las particulas en la poblacion*/
double **pbestsFit;	 /* valor de aptitud pBest de las particulas en la poblacion*/
double **velocity;	 /* velocidad de las particulas en la poblacion*/
double pMut = 0.5;	 /* probabilidad de mutacion*/
unsigned int nondomCtr = 0;  /* numero de soluciones no dominadas en el archivo*/

/*
Variables para el calculo del hipervolumen
*/
double *contribution;
double *points;
double *contribution_2;
double *contribution_3;

double *pContribution;
double *bpContribution;
//double *used;
double *referencePoints;
double *referencePoints_2;


double delta = 10.0;
unsigned int param_k;
double lowerbound = 0.0;
double upperbound;
double volumeT;
double *volumeFit;
double **archiveFit_2;
//double **archiveFit_3;
FILE *hv;

void initialize_memory(void);
void free_memory(void);
void initialize_data(unsigned int,unsigned int);
void initialize_rand(void);
void initialize_pop(void);
void initialize_vel(void);
void evaluate(void);

/*Test function*/
double kita_f1(unsigned int);
double kita_f2(unsigned int);
double kursawe_f1(unsigned int);
double kursawe_f2(unsigned int);
double deb_f1(unsigned int);
double deb_f2(unsigned int);

void deb2(unsigned int);
void deb3(unsigned int);
void fonseca2(unsigned int);

double DTLZ6_f1(unsigned int);
double DTLZ6_f2(unsigned int);
double DTLZ6_f3(unsigned int);

void ZDT1(unsigned int);
void ZDT2(unsigned int);
void ZDT3(unsigned int);
void ZDT4(unsigned int);
void ZDT6(unsigned int);

void DTLZ1(unsigned int);
void DTLZ2(unsigned int);
void DTLZ3(unsigned int);
void DTLZ4(unsigned int);
void DTLZ5(unsigned int);
void DTLZ6(unsigned int);
void DTLZ7(unsigned int);

void store_pbests(void);
void insert_nondom(void);
void delete_particle(unsigned int);
unsigned int compare_particles(unsigned int, unsigned int);
unsigned int check_constraints(double *);
void mutate(unsigned int);
void get_ranges(double *, double *);
void maintain_particles(void);
unsigned int compute_gBest(void);
void compute_hvContribution(void);

void compute_velocity(void);
void compute_position(void);
void update_archive(void);
unsigned int check_nondom(unsigned int);
void update_pbests(void);
void save_results(char *);


/*Algoritmo Principal*/
int main(int argc, char *argv[])
{
    char name[20];// = "_kita_1";
    char archiveName[20] = "archive";
	char fitputName[20] = "fitput";
	char varputName[20] = "varput";
	char hvName[20]= "hv";
	unsigned int i, j, t;
	unsigned int fun, gen;
	clock_t  startTime, endTime;
	double duration, clocktime;
	FILE *fitfile,*varfile;
/*Parametros iniciales*/

    strcpy(name,argv[1]);
    fun = atoi(argv[2]);
    gen = atoi(argv[3]);

    printf("%s %d %d \n",name,fun,gen);

    initialize_data(fun,gen);
/*Iniciar variables*/
	initialize_memory();

	/*
	sprintf(name,"kita_1_");
	sprintf(archiveName,strcat(name,"archive.out"));
	sprintf(fitputName, strcat(name,"fitput.out"));
	sprintf(varputName, strcat(name,"varput.out"));
	sprintf(hvName,strcat(name,"hv.out"));
*/
	//fitfile = fopen(strcat(strcat(fitputName,name),".out"),"w");
	//varfile = fopen(strcat(strcat(varputName,name),".out"),"w");
	//hv = fopen(strcat(strcat(hvName,name),".out"),"w");
/* Iniciar generador de aleatorios*/
	initialize_rand();
	startTime = clock();
/* Iniciar contador de generaciones*/
	t = 0;
/* Iniciar valores de la poblacion de manera aleatoria*/
	initialize_pop();
/* Calcular velocidad inicial*/
	initialize_vel();
/* Evaluar las particulas de la poblacion*/
	evaluate();
/* Almacenar el pBest inicial (variables and valor de aptitud) de las particulas*/
	store_pbests();
/* Insertar las particulas no domindas de la poblacion en el archivo*/

	insert_nondom();
	//printf("\n%d",nondomCtr);
/*Ciclo Principal*/
	while(t <= maxgen)
	{
		clocktime = (clock() - startTime)/(double)CLOCKS_PER_SEC;
		if(verbose > 0 && t%printevery==0 || t == maxgen)
		{
			fprintf(stdout,"%d  %.2f \n",t,clocktime);
			fflush(stdout);
		}
		/*if(t%printevery==0 || t == maxgen)
		{
			fprintf(outfile,"Generation %d Time: %.2f sec\n",t,clocktime);
		}*/
    /* Calcular la nueva velocidad de cada particula en la pooblacion*/
		//printf("\n 1");
		compute_velocity();

    /* Calcular la nueva posicion de cada particula en la poblacion*/
        //printf("\n 2");
        compute_position();
    /* Mantener las particulas en la poblacion de la poblacion en el espacio de busqueda*/
		//printf("\n 3");
		maintain_particles();
    /* Pertubar las particulas en la poblacion*/
		if(t < maxgen * pMut)
			mutate(t);

    /* Evaluar las particulas en la poblacion*/
		//printf("\n 4");
		evaluate();
    /* Insertar nuevas particulas no domindas en el archivo*/
		//printf("\n 5");
		update_archive();
    /* Actualizar el pBest de las particulas en la poblacion*/
		//printf("\n 6");
		update_pbests();
    /* Escribir resultados del mejor hasta ahora*/
    //verbose > 0 && t%printevery==0 || t == maxgen
		/*if(t%printevery==0 || t == maxgen)
		{
			//fprintf(outfile, "Size of Pareto Set: %d\n", nondomCtr);
			fprintf(fitfile, "%d\n",t);
			fprintf(varfile, "%d\n",t);
			for(i = 0; i < nondomCtr; i++)
			{
			    for(j = 0; j < maxfun; j++)
                    fprintf(fitfile, "%f ", archiveFit[i][j]);
                fprintf(fitfile, "\n");
			}
			fprintf(fitfile, "\n\n");
			fflush(fitfile);
			for(i = 0; i < nondomCtr; i++)
			{
			    for(j = 0; j < maxfun; j++)
                    fprintf(varfile, "%f ", archiveVar[i][j]);
                fprintf(varfile, "\n");
			}
			fprintf(varfile, "\n\n");
			fflush(varfile);
		}*/
    /* Incrementar contador de generaciones*/
        t++;
        //printf("%d\n",t);
	}
 /* Escribir resultados en el archivo */
	save_results(strcat(strcat(archiveName,name),".out"));
	endTime = clock();
	duration = ( endTime - startTime ) / (double)CLOCKS_PER_SEC;
	fprintf(stdout, "%lf sec\n", duration);
	//fclose(fitfile);
	//fclose(varfile);
	free_memory();
	return EXIT_SUCCESS;
}
void initialize_data(unsigned int fun,unsigned int mg)
{
    function = fun;
    maxgen = mg;         /* numero maximo de generaciones */

    switch(fun)
	{
		case 100: /* Kita */
            maxfun = 2;
            maxvar = 2;
            optimization = 1;
            break;
		case 200: /* Kursawe */
			maxfun = 2;
            maxvar = 3;
            optimization = 0;
			break;
		case 300: /* Deb 1 */
			maxfun = 2;
            maxvar = 2;
            optimization = 0;
			break;
        case 350: /* Deb 2 */
			maxfun = 2;
            maxvar = 2;
            optimization = 0;
			break;
        case 400: /* Deb 3 */
			maxfun = 2;
            maxvar = 2;
            optimization = 0;
			break;
        case 450: /* Fonseca 2 */
			maxfun = 2;
            maxvar = 2;
            optimization = 0;
			break;
        case 500: /* DTLZ7 */
			maxfun = 3;
            maxvar = 22;
            optimization = 0;
            break;
        case 600: /* zdt1 */
			maxfun = 2;
            maxvar = 30;
            optimization = 0;
            break;
        case 605: /* zdt2 */
			maxfun = 2;
            maxvar = 30;
            optimization = 0;
            break;
        case 610: /* zdt3 */
			maxfun = 2;
            maxvar = 30;
            optimization = 0;
            break;
        case 615: /* zdt4 */
			maxfun = 2;
            maxvar = 10;
            optimization = 0;
            break;
        case 620: /* zdt6 */
			maxfun = 2;
            maxvar = 10;
            optimization = 0;
            break;
        case 700: /* DTLZ1*/
			maxfun = 3;
            maxvar = 12;
            optimization = 0;
            break;
        case 705: /* DTLZ2*/
			maxfun = 3;
            maxvar = 12 ;
            optimization = 0;
            break;
        case 710: /* DTLZ3*/
			maxfun = 3;
            maxvar = 12;
            optimization = 0;
            break;
        case 715: /* DTLZ4*/
			maxfun = 3;
            maxvar = 12;
            optimization = 0;
            break;
        case 720: /* DTLZ5*/
			maxfun = 3;
            maxvar = 12;
            optimization = 0;
            break;
        case 725: /* DTLZ6*/
			maxfun = 3;
            maxvar = 12;
            optimization = 0;
            break;
        case 730: /* DTLZ7*/
			maxfun = 3;
            maxvar = 22;
            optimization = 0;
            break;
             /** Agregar mas aqui ***/
	}
}

void initialize_memory()
{
	unsigned int i, j;
	/*reservando Memoria*/
	archiveVar = (double **) calloc(archive_size,sizeof(double *));
	archiveFit = (double **) calloc(archive_size,sizeof(double *));
	archiveFit_2 = (double **) calloc(archive_size,sizeof(double *));
	//archiveFit_3 = (double **) calloc(archive_size,sizeof(double *));
	for(i=0;i < archive_size;i++)
	{
		archiveVar[i] = (double *) calloc(maxvar,sizeof(double));
		archiveFit[i] = (double *) calloc(maxfun,sizeof(double));
		archiveFit_2[i] = (double *) calloc(maxfun,sizeof(double));
		//archiveFit_3[i] = (double *) calloc(maxfun,sizeof(double));
	}

	popVar = (double **) calloc(popsize,sizeof(double *));
	popFit = (double **) calloc(popsize,sizeof(double *));
	pbestsVar = (double **) calloc(popsize,sizeof(double *));
	pbestsFit = (double **) calloc(popsize,sizeof(double *));
	velocity = (double **) calloc(popsize,sizeof(double *));

	for(i=0;i < popsize;i++)
	{
		popVar[i] = (double *) calloc(maxvar,sizeof(double));
		pbestsVar[i] = (double *) calloc(maxvar,sizeof(double));
		velocity[i] = (double *) calloc(maxvar,sizeof(double));
	}
	for(i=0;i < popsize;i++)
	{
		popFit[i] = (double *) calloc(maxfun,sizeof(double));
		pbestsFit[i] = (double *) calloc(maxfun,sizeof(double));
	}
	referencePoints = (double *) calloc(maxfun,sizeof(double));
	referencePoints_2 = (double *) calloc(maxfun,sizeof(double));
	//volumeFit = (double *)calloc(archive_size,sizeof(double));

	contribution = (double *)calloc(archive_size,sizeof(double));
    points = (double *)calloc(archive_size * maxfun,sizeof(double));

    contribution_2 = (double *)calloc(archive_size,sizeof(double));
    contribution_3 = (double *)calloc(archive_size,sizeof(double));

    pContribution = (double *)calloc(archive_size,sizeof(double));
    bpContribution = (double *)calloc(archive_size,sizeof(double));
}

void free_memory()
{
    int i;
	for(i=0;i < archive_size;i++)
	{
		free(archiveVar[i]);
		free(archiveFit[i]);
		free(archiveFit_2[i]);
	//	free(archiveFit_3[i]);
	}
    free(archiveVar);
	free(archiveFit);
	free(archiveFit_2);
	//free(archiveFit_3);

    for(i=0;i < popsize;i++)
	{
		free(popVar[i]);
		free(pbestsVar[i]);
		free(velocity[i]);
	}
	for(i=0;i < popsize;i++)
	{
		free(popFit[i]);
		free(pbestsFit[i]);
	}
	free(popVar);
	free(popFit);
	free(pbestsVar);
	free(pbestsFit);
	free(velocity);

	free(referencePoints);
	free(referencePoints_2);
	free(volumeFit);

	free(points);
	free(contribution);
	free(contribution_2);
	free(contribution_3);

	free(pContribution);
    free(bpContribution);
}

void initialize_rand() /* Iniciando el generador de aleatorios randomlib.h */
{
	unsigned int i, j;
	srand((unsigned int)time((time_t *)NULL));
	i = (unsigned int) (31329.0 * rand() / (RAND_MAX + 1.0));
	j = (unsigned int) (30082.0 * rand() / (RAND_MAX + 1.0));
	RandomInitialise(i,j);
}

void initialize_pop() /* Iniciando la poblacion */
{
	unsigned int i, j;
	switch(function)
	{
		case 100: /* Kita */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 7.0);
			break;
		case 200: /* Kursawe */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(-5.0, 5.0);
			break;
		case 300: /* Deb 1 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.1, 0.8191);
			break;
        case 350: /* Deb 2 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);// 0.8191);
			break;
        case 400: /* Deb 3 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
        case 450: /* Fonseca 2 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(-4.0, 4.0);
			break;

		case 500: /* DTLZ6 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
        case 600: /* zdt1 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(-0.0, 1.0);
			break;
        case 605: /* zdt2 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
        case 610: /* zdt3 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
        case 615: /* zdt4 */
			for(i=0; i < popsize; i++)
			{
			    popVar[i][0] = RandomDouble(0.0, 1.0);
			    for(j =1; j < maxvar; j++)
					popVar[i][j] = RandomDouble(-5.0, 5.0);
			}
			break;
        case 620: /* zdt6 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
        case 700: /* dtlz1 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
        case 705: /* dtlz2 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
        case 710: /* dtlz3 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
        case 715: /* dtlz4 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
        case 720: /* dtlz5 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
        case 725: /* dtlz6 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
        case 730: /* dtlz7 */
			for(i=0; i < popsize; i++)
				for(j =0; j < maxvar; j++)
					popVar[i][j] = RandomDouble(0.0, 1.0);
			break;
    /** Agregar maa aqui ***/
	}
}

void initialize_vel() /* Velocidad inicial igual a cero */
{
	unsigned int i, j;
	for(i = 0; i < popsize; i++)
		for(j = 0; j < maxvar; j++)
			velocity[i][j] = 0.0;
}

void evaluate() /* Evaluar las particulas de la poblacion */
{
	unsigned int i, j;
	for(i = 0; i < popsize; i++)
	{
		for(j = 0; j < maxfun; j++)
		{
		    //
			switch(function + j)
			{

				case 100:	/* Kita's primer objetivo */
					popFit[i][j] = kita_f1(i);
					break;
				case 101:	/* Kita's segundo objetivo */
					popFit[i][j] = kita_f2(i);
					break;
				case 200:	/* Kursawe's primer objetivo */
					popFit[i][j] = kursawe_f1(i);
					break;
				case 201:	/* Kursawe's segundo objetivo */
					popFit[i][j] = kursawe_f2(i);
					break;
				case 300:	/* Deb's primer objetivo */
					popFit[i][j] = deb_f1(i);
					break;
				case 301:	/* Deb's segundo objetivo */
					popFit[i][j] = deb_f2(i);
					break;
				case 350:
				//printf("%d %d\n",i,function);
                    deb2(i);
                    break;
				case 400:
                    deb3(i);
                    break;
                case 450:
                    fonseca2(i);
                    break;
				case 500:	/* DTLZ6's primer objetivo */
					popFit[i][j] = DTLZ6_f1(i);
					break;
				case 501:	/* DTLZ6's segundo objetivo */
					popFit[i][j] = DTLZ6_f2(i);
					break;
				case 502:	/* DTLZ1's tercer objetivo */
					popFit[i][j] = DTLZ6_f3(i);
					break;
                case 600:	/* zdt1 primer objetivo*/
					ZDT1(i);
					break;
                case 605:	/* zdt2 primer objetivo*/
					ZDT2(i);
					break;
                case 610:	/* zdt3 primer objetivo*/
					ZDT3(i);
					break;
                case 615:	/* zdt4 primer objetivo*/
                    ZDT4(i);
					break;
                case 620:	/* zdt6 segundo objetivo*/
					ZDT6(i);
					break;
                case 700:
                    DTLZ1(i);
                    break;
                case 705:
                    DTLZ2(i);
                    break;
                case 710:
                    DTLZ3(i);
                    break;
                case 715:
                    DTLZ4(i);
                    break;
                case 720:
                    DTLZ5(i);
                    break;
                case 725:
                    DTLZ6(i);
                    break;
                case 730:
                    DTLZ7(i);
                    break;
	/** Agregar mas aqui **/
			} /* end of switch */
		}
	}
}

void store_pbests() /* Almacenar pBest (variables y valores de fitness) */
{
	unsigned int i, j;
  /* Alamcenar los valores de las variable */
	for(i=0; i < popsize; i++)
		for (j = 0; j < maxvar; j++)
			pbestsVar[i][j] = popVar[i][j];
  /* Almacenar el fitness */
	for(i=0; i < popsize; i++)
		for (j = 0; j < maxfun; j++)
			pbestsFit[i][j] = popFit[i][j];
	for (j = 0; j < maxfun; j++)
        referencePoints_2[j]=RMIN;
}

void insert_nondom() /* Insertar particulas no dominadas en el archivo */
{
	unsigned int i, j, k, total, insertFlag, bottom, cont;
	double archiveCons[maxvar], popCons[maxvar];

	for(i=0; i< popsize; i++)
	{
		if(nondomCtr == 0)
		{ /* Si el archivo esta vacio */
      /* Insertar la particula en el archivo */
			for(j = 0; j < maxvar; j++)
				archiveVar[nondomCtr][j] = popVar[i][j];
			for(j = 0; j < maxfun; j++)
				archiveFit[nondomCtr][j] = popFit[i][j];
			nondomCtr += 1;
		}
		else
		{ /* Si el archivo no esta vacio */
			insertFlag = 1;
            /*Para cada particula en el archivo */
			for(k=0; k < nondomCtr; k++)
			{
            /* Primero, verificar si es factibe */
				for(j=0; j < maxvar; j++)
				{
					popCons[j] = popVar[i][j];
					archiveCons[j] = archiveVar[k][j];
				}
                    /* Si ambas partculas no son factibles */
				if((check_constraints(archiveCons) > 0) && (check_constraints(popCons) > 0))
				{
					delete_particle(k); /* Eliminar particula del archivo */
					insertFlag = 0;		/* No insertar la particula en la poblacion */
					break;
				}
				else if(check_constraints(popCons) > 0)
				{ /* Si la particula de la poblacion en infactible */
					insertFlag = 0;/* No insertar la particula en la poblacion */
					break;
				}
				else if(check_constraints(archiveCons) > 0)
				{ /* Si la particula en el archivo es infactible */
					delete_particle(k); /* Borrar la particula en el archivo  */
					if((nondomCtr != 0) || (k != nondomCtr-1))
						k--;
					continue;
				}
	/* Segundo, verificar dominancia*/
				total = 0;
				cont = 0;
                /* Si ambos son factiblesf, vereficar no diminancia*/
				for(j=0; j < maxfun; j++)
				{
					if(( (popFit[i][j] < archiveFit[k][j]) && (optimization == 0)) || (popFit[i][j] > archiveFit[k][j]) && (optimization == 1))
                        total += 1;
				}
				for(j=0; j < maxfun; j++)
				{
					if( (popFit[i][j] == archiveFit[k][j]))
                    {
                        cont ++;     /* No insertar la particula en la poblacion */
                        insertFlag = 0;
                    }
				}
				if(total == maxfun || cont == maxfun)   /* Si la particula es totalmente dominada     */
                    delete_particle(k); /* Borrar particula del archivo*/
				else if(total == 0)
				{  /* Si la particula en dominada en el archivo*/
					insertFlag = 0;     /* No insertar la particula en la poblacion */
					break;
				}

			} /* Termina la comparacion de una particula en la poblacion con el archivo */
		}
    /* Insertar la particula si es factible y no dominada */
		if(insertFlag == 1)
		{
            /* Si la memoria no esta llena, insertar particula */
      		if (nondomCtr < archive_size)
			{
				for(j = 0; j < maxvar; j++)
					archiveVar[nondomCtr][j] = popVar[i][j];
				for(j = 0; j < maxfun; j++)
					archiveFit[nondomCtr][j] = popFit[i][j];
				nondomCtr += 1;
			}
			else
			{   /* Calcular hipervolumen*/

				bottom = (unsigned int)((nondomCtr-1) * 0.90);
            /* Selecionar aleatoriamente un lugar para reemplazar */
				k = RandomInt(bottom, nondomCtr-1);
            /* Insertar nueva particula en el archivo */
				for(j = 0; j < maxvar; j++)
					archiveVar[k][j] = popVar[i][j];
				for(j = 0; j < maxfun; j++)
					archiveFit[k][j] = popFit[i][j];
			}
		}
	} /* Termina la comparacion de particulas de la poblacion en el archivo */
}

void delete_particle(unsigned int k) /* Eliminar Particula en el archivo */
{
	unsigned int j;

	/* Si la particula es infactible o la ultima particula del archivo o solo una particula en el archivo */
	if( (nondomCtr == 1) || (k == (nondomCtr-1)) )
	{
		nondomCtr -= 1;
	}
	else
	{       /* mover la ultima particula en el lugar de la particula infactible */
		for(j = 0; j < maxvar; j++)
			archiveVar[k][j] = archiveVar[nondomCtr-1][j];
		for(j = 0; j < maxfun; j++)
			archiveFit[k][j] = archiveFit[nondomCtr-1][j];
		nondomCtr -= 1;
	}
}

unsigned int check_constraints(double *consVar) /* Verificar restricciones q detereminan la factibilidad*/
{
	unsigned int violations = 0;
	switch(function)
	{
		case 100: /* Kita  */
			if( 0 < (double)(consVar[0] / 6.0 + consVar[1] - 13.0 / 2.0 ))
				violations++;
			if( 0 < (double)(consVar[0] / 2.0 + consVar[1] - 15.0 / 2.0))
				violations++;
			if( 0 < (double)(5.0 * consVar[0] + consVar[1] - 30.0) )
				violations++;
			return violations;
    /** Agregar mas aqui**/
		default:
			return 0;
	}
}

void compute_hvContribution()
{
    /*Calcula al gBest conforme a la mejor contribucion*/
    int i, j, k;
    double best;
    //contribution = (double *)calloc(nondomCtr,sizeof(double));
    //points = (double *)calloc(nondomCtr * maxfun,sizeof(double));
    param_k = 1;

    /*for(i=0; i < nondomCtr; i++)
      {
          for(j=0; j < maxfun; j++)
        {
            printf("%f ",archiveFit[i][j]);
        }
         printf("\n");
      }*/
    /*printf("\n");
    for(i=0, k = 0; i < nondomCtr; i++)
        for(j=0; j < maxfun; j++)
        {
            points[k] = archiveFit[i][j];
            k++;
        }
    printf("\n %d",nondomCtr);
    //printf("\n Puntos");
    for(i=0;i< (nondomCtr*maxfun);i++)
    {
        printf("\n %f",points[i]);
    }*/

    for(j=0; j < maxfun; j++)
    {
        referencePoints[j] = popFit[0][j];
        for(i=0; i < nondomCtr; i++)
        {
            if(popFit[i][j] > referencePoints[j])
                 referencePoints[j] = popFit[i][j];
        }
    }
    for(j=0; j < maxfun; j++)
    {
        if((referencePoints[j] + delta) > referencePoints_2[j])
        {
            referencePoints[j] = referencePoints[j] + delta;
            referencePoints_2[j] = referencePoints[j] + delta;
        }
        else
        {
             referencePoints[j]=referencePoints_2[j];
        }

    }

   /* for(j=0; j < maxfun; j++)
    {
        referencePoints[j] = bound;
    }*/
    if(maxfun > 2)
    {
        upperbound = referencePoints[0];
        for(j=0; j < maxfun; j++)
            if(referencePoints[j] > upperbound)
                upperbound = referencePoints[j];

        referencePoints[0] = lowerbound;
        //referencePoints[1] = upperbound + delta;
        referencePoints[1] = upperbound;
    }
    hypeIndicator(contribution,nondomCtr,referencePoints,param_k,points,maxfun);

    /*for(i=0;i< nondomCtr;i++)
    {
        printf("\n %f",contribution[i]);
    }*/

}

unsigned int compute_gBest()
{
    /*Como calcular el gBest para cada particula */
    unsigned int k, i,j;
    double best;
    best = contribution[0];
    k = 0;
    for(i=1; i < nondomCtr; i++)
    {
        if(contribution[i] > best)
        {
            best = contribution[i];
            k = i;
        }

    }
    //printf("\n %e %d",best, k);

   // free(contribution);
   // free(points);
    return k;
}

void compute_velocity() /* Calcular la nueva velocidad de cada particula en la poblacion */
{
	unsigned int i, j, gBest;
    compute_hvContribution();
    gBest=compute_gBest();
	for(i = 0; i < popsize; i++)
	{
		for(j = 0; j < maxvar; j++)
			velocity[i][j] = 0.4 * velocity[i][j] + 1.0 * RandomDouble(0.0, 1.0) * (pbestsVar[i][j] - popVar[i][j]) + 1.0 * RandomDouble(0.0, 1.0) * (archiveVar[gBest][j] - popVar[i][j]);
            //velocity[i][j] = 1.0 * velocity[i][j] + 1.0 * RandomDouble(0.0, 1.0) * (pbestsVar[i][j] - popVar[i][j]) + 1.0 * RandomDouble(0.0, 1.0) * (archiveVar[gBest][j] - popVar[i][j]);
            /* W* Vi + C1  * RandomDouble(0.0, 1.0) * (pBest - Xi) + C2  * RandomDouble(0.0, 1.0) * (gBest  -  Xi) */
	}
}

void compute_position()
{
    unsigned int i, j;
    /* Calcular la nueva posicion de las particulas */
	for(i = 0; i < popsize; i++)
		for(j = 0; j < maxvar; j++)
            popVar[i][j] = popVar[i][j] + velocity[i][j];
}

void mutate(unsigned int t) /* Mutacion de PArticulas MOPSO */
{
	unsigned int i;
	int dimension = 0;
	double minvalue[maxvar], maxvalue[maxvar], minvaluetemp, maxvaluetemp, range;
	double valtemp = 0;
	get_ranges(minvalue,maxvalue);
	for(i=0; i < popsize; i++)
	{
		if(flip(pow(1.0-(double)t / (maxgen * pMut),1.5)))
		{
			dimension = RandomInt(0,maxvar-1);
			range =(maxvalue[dimension] - minvalue[dimension]) * pow(1.0 - (double)t / (maxgen * pMut),1.5) / 2;
			valtemp = RandomDouble(range, -range);
			if(popVar[i][dimension] - range < minvalue[dimension])
				minvaluetemp = minvalue[dimension];
			else
				minvaluetemp = popVar[i][dimension] - range;
			if(popVar[i][dimension] + range > maxvalue[dimension])
				maxvaluetemp = maxvalue[dimension];
			else
				maxvaluetemp = popVar[i][dimension] + range;
			popVar[i][dimension] = RandomDouble(minvaluetemp, maxvaluetemp);
		}
	}
}

void update_archive()  /* Insertar nuevas particulas no dominadas de la poblacion en el archivo */
{
	unsigned int i, j, k, bottom, m, l, n, flag, temp, pWorse, cont;
    double tempE, worse;
    if(acotado == 1)
    {
        for(i = 0, n = 0; i < nondomCtr; i++)
            for(j = 0; j < maxfun; j++)
            {
                points[n] = archiveFit[i][j];
                archiveFit_2[i][j] = archiveFit[i][j];
                n++;
            }
        hypeIndicator(contribution_2,nondomCtr,referencePoints,param_k,points,maxfun);
        worse= contribution_2[0];
        pWorse = 0;
        for(i = 1; i < nondomCtr; i++)
        {
            if(contribution_2[i] < worse)
            {
                worse = contribution_2[i];
                pWorse = i;
            }
        }
    }
  /* para cada particula en la poblacion */
	for(k = 0; k < popsize; k++)
	{
    /* Si la particula en la poblacion es no dominada */
       // cont = 0;
		if (check_nondom(k) == 1)
		{
            /* Si la memoria no esta llena, insertar particula */
			if (nondomCtr < archive_size)
			{
			    i = nondomCtr;

				for(j = 0; j < maxvar; j++)
					archiveVar[i][j] = popVar[k][j];
				for(j = 0; j < maxfun; j++)
					archiveFit[i][j] = popFit[k][j];
				nondomCtr += 1;
            }
			else
			{      // Calcular hipervolumen
			    if (acotado == 1)
                {
                    for(j = 0; j < maxfun; j++)
                        archiveFit_2[pWorse][j] = popFit[k][j];
                    for(i = 0, n = 0; i < nondomCtr; i++)
                        for(j = 0; j < maxfun; j++)
                        {
                            points[n] = archiveFit_2[i][j];
                            n++;
                        }
                    hypeIndicator(contribution_2,nondomCtr,referencePoints,param_k,points,maxfun);
                    if(contribution_2[pWorse] > worse)
                    {
                        for(j = 0; j < maxvar; j++)
                            archiveVar[pWorse][j] = popVar[k][j];
                        for(j = 0; j < maxfun; j++)
                            archiveFit[pWorse][j] = popFit[k][j];
                    }
                    for(i = 0, n = 0; i < nondomCtr; i++)
                        for(j = 0; j < maxfun; j++)
                        {
                            points[n] = archiveFit[i][j];
                            n++;
                        }
                    hypeIndicator(contribution_2,nondomCtr,referencePoints,param_k,points,maxfun);
                    worse= contribution_2[0];
                    pWorse = 0;
                    for(i = 1; i < nondomCtr; i++)
                    {
                        if(contribution_2[i] < worse)
                        {
                            worse = contribution_2[i];
                            pWorse = i;
                        }
                    }
                    /*for(j = 0; j < maxvar; j++)
                        archiveVar[pWorse][j] = popVar[k][j];
                    for(j = 0; j < maxfun; j++)
                        archiveFit[pWorse][j] = popFit[k][j];
                    for(i = 0, n = 0; i < nondomCtr; i++)
                        for(j = 0; j < maxfun; j++)
                        {
                            points[n] = archiveFit[i][j];
                            n++;
                        }
                    hypeIndicator(contribution_2,nondomCtr,referencePoints,param_k,points,maxfun);
                    worse= contribution_2[0];
                    pWorse = 0;
                    for(i = 1; i < nondomCtr; i++)
                    {
                        if(contribution_2[i] < worse)
                        {
                            worse = contribution_2[i];
                            pWorse = i;
                        }
                    }*/
                }
                else
                {
                    bottom = (unsigned int)((nondomCtr-1) * 0.90);
            /* Selecionar aleatoriamente un lugar para reemplazar */
                    l = RandomInt(bottom, nondomCtr-1);
                    for(j = 0; j < maxvar; j++)
                        archiveVar[l][j] = popVar[k][j];
                    for(j = 0; j < maxfun; j++)
                        archiveFit[l][j] = popFit[k][j];
                }
			}
		}
		//printf("calulando extremos\n");
	}
	if(function == 620)
        for(i = 0; i < nondomCtr; i++)
        {
            if(archiveVar[i][1] > 1.0)
            {       /* mover la ultima particula en el lugar de la particula infactible */
                for(j = 0; j < maxvar; j++)
                    archiveVar[i][j] = archiveVar[nondomCtr-1][j];
                for(j = 0; j < maxfun; j++)
                    archiveFit[i][j] = archiveFit[nondomCtr-1][j];
                nondomCtr -= 1;
            }
        }
	//printf("aaa\n");
}

unsigned int check_nondom(unsigned int i) /* Check for feasibility and nondomination of particles in pop and archive */
{
	unsigned int sum, h = 0, j;//, cont;
	double archiveCons[maxvar], popCons[maxvar];
	do
	{
		sum = 0;
		//cont = 0;
		for(j=0; j < maxvar; j++)
		{
			archiveCons[j] = archiveVar[h][j];
			popCons[j] = popVar[i][j];
		}
        /* Si la particula en el archivo tiene menos restricciones violadas */
		if(check_constraints(archiveCons) < check_constraints(popCons))
			return 0;
        /* Si la particula tiene mas restriccion, eliminar*/
		if(check_constraints(archiveCons) > check_constraints(popCons))
		{
			for(j = 0; j < maxvar; j++)
				archiveVar[h][j] = archiveVar[nondomCtr-1][j];
			for(j = 0; j < maxfun; j++)
				archiveFit[h][j] = archiveFit[nondomCtr-1][j];
			nondomCtr -= 1;
		}
		else /*Si tienen las mismas violacion de restricciones*/
		{
			for(j = 0; j < maxfun; j++)
			{
				if( ((archiveFit[h][j] < popFit[i][j]) && (optimization == 0)) || ((archiveFit[h][j] > popFit[i][j]) && (optimization == 1)))
                    sum += 1;
              //  if( archiveFit[h][j] == popFit[i][j])
                //    cont += 1;
			}
			if(sum == maxfun)// || cont == maxfun)	/* Si la particula en el archio es dominada */
				return 0;
			else if(sum == 0)
			{	/* Si la particula en el archivo es domianda borrar *///Revisar
				for(j = 0; j < maxvar; j++)
					archiveVar[h][j] = archiveVar[nondomCtr-1][j];
				for(j = 0; j < maxfun; j++)
					archiveFit[h][j] = archiveFit[nondomCtr-1][j];
				nondomCtr -= 1;
			}
			else
			{
				h += 1;
			}
		}
	} while(h < nondomCtr);
	return 1;
}

void get_ranges(double *minvalue,double *maxvalue) /* Obtener rango de valores de las variables */
{
	unsigned int i;
	switch (function)
	{
		case 100:	/* Kita */
		//printf("Arreglando Intervalos de las variable\n");
		for(i = 0; i < maxvar; i++)
		{
			minvalue[i] = 0.0;
			maxvalue[i] = 7.0;
		}
		break;
	case 200:   /* Kursawe */
		for(i = 0; i < maxvar; i++)
		{
			minvalue[i] = -5.0;
			maxvalue[i] = 5.0;
		}
		break;
	case 300:
      	for(i = 0; i < maxvar; i++)
		{
			minvalue[i] = 0.1;
			maxvalue[i] = 0.8191;
		}
		break;

	case 350: case 400:	/* Deb 2*/
        //printf("Arreglando Intervalos de las variable\n");
		for(i = 0; i < maxvar; i++)
		{
			minvalue[i] = 0.0;
			maxvalue[i] = 1.0;
		}
		break;
	case 450:	/* Fonseca 2 */
		for(i = 0; i < maxvar; i++)
		{
			minvalue[i] = -4.0;
			maxvalue[i] = 4.0;
		}
		break;
	case 500:	/* DTLZ6 */
		for(i = 0; i < maxvar; i++)
		{
			minvalue[i] = 0.0;
			maxvalue[i] = 1.0;
		}
		break;
    case 600: case 605: case 610: /* zdt*/
		for(i = 0; i < maxvar; i++)
		{
			minvalue[i] = -0.0;
			maxvalue[i] = 1.0;
		}
		break;
    case 615: /*zdt4*/
        minvalue[0] = 0.0;
        maxvalue[0] = 1.0;
		for(i=1; i < maxvar; i++)
			{
			    minvalue[i] = -5.0;
			    maxvalue[i] =  5.0;
			}
        break;
    case 620:	/* zdt6*/
		for(i = 0; i < maxvar; i++)
		{
			minvalue[i] = 0.0;
			maxvalue[i] = 1.0;
		}
		break;
    case 700: case 705: case 710: case 715: case 720:	/* dtlz*/
    case 725: case 730:
		for(i = 0; i < maxvar; i++)
		{
			minvalue[i] = 0.0;
			maxvalue[i] = 1.0;
		}
		break;

	}
}

void maintain_particles()  /* Mantener las particulas de la poblacion en el espacio de busqueda */
{
	unsigned int i, j;
	double minvalue[maxvar], maxvalue[maxvar];
	//return;
    /* Verifica la estabilidad de la poblacion al mantener las partículas en sus fronteras */
    switch (function)
	{
		case 100: /* Kita */
            for(i = 0; i < maxvar; i++)
			{
				minvalue[i] = 0.0;
				maxvalue[i] = 7.0;
			}
			break;
		case 200: /* Kursawe */
			for(i = 0; i < maxvar; i++)
			{
				minvalue[i] = -5.0;
				maxvalue[i] = 5.0;
			}
			break;
		case 300:
            for(i = 0; i < maxvar; i++)
			{
				minvalue[i] = 0.1;
				maxvalue[i] = 0.8191;
			}
			break;
		case 350: case 400:/* Deb */
			for(i = 0; i < maxvar; i++)
			{
				minvalue[i] = 0.0;
				maxvalue[i] = 1.0;
			}
			break;
        case 450:	/* Fonseca 2 */
            for(i = 0; i < maxvar; i++)
            {
                minvalue[i] = -4.0;
                maxvalue[i] = 4.0;
            }
		break;
		case 500: /* DTLZ6 */
			for(i = 0; i < maxvar; i++)
			{
				minvalue[i] = 0.0;
				maxvalue[i] = 1.0;
			}
			break;
        case 600: case 605: case 610: /* zdt*/
            for(i = 0; i < maxvar; i++)
            {
                minvalue[i] = -0.0;
                maxvalue[i] = 1.0;
            }
            break;
        case 615:
            minvalue[0] = 0.0;
            maxvalue[0] = 1.0;
            for(i=1; i < maxvar; i++)
                {
                    minvalue[i] = -5.0;
                    maxvalue[i] =  5.0;
                }
            break;
        case 620:	/* zdt6*/
            for(i = 0; i < maxvar; i++)
            {
                minvalue[i] = 0.0;
                maxvalue[i] = 1.0;
            }
            break;
        case 700: case 705: case 710: case 715: case 720:	/* zdt*/
        case 725: case 730:
            for(i = 0; i < maxvar; i++)
            {
                minvalue[i] = 0.0;
                maxvalue[i] = 1.0;
            }
            break;

  /** Agregar mas aqui! **/
	}
	for(i = 0; i < popsize; i++)
	{
		for(j = 0; j < maxvar; j++)
		{
            /* Si la particula se va mas lejos del valor minimo */
			if(popVar[i][j] < minvalue[j])
			{
                /* HAcer igual al valor minimo */
				popVar[i][j] = minvalue[j];
                /* Cambiar la direccion de manera opuesta */
				velocity[i][j] = -velocity[i][j];
				/*cambiar la velocidad de la particula a cero*/
				//velocity[i][j] = 0;
            }
            /* Si la particula se va mas lejos del valor maximo */
            if(popVar[i][j] > maxvalue[j])
            {
                /* Hacer igual al valor maximo*/
                popVar[i][j] = maxvalue[j];
                /* Cambiar la direccion de manera opuesta */
				velocity[i][j] = -velocity[i][j];
				/*cambiar la velocidad de la particula a cero*/
				//velocity[i][j] = 0;
            }
        }
    }
}

void update_pbests() /* Update personal bests of particles in the population */
{
	unsigned int i,n, j, sum, better;

	/*for(i = 0, n= 0; i < nondomCtr; i++)
        for(j = 0; j < maxfun; j++)
        {
            points[n] = popFit[i][j];
            n++;
        }
    hypeIndicator(pContribution,popsize,referencePoints,param_k,points,maxfun);*/

	for(i = 0; i < popsize; i++)
        {
            sum = 0;
            for(j = 0; j < maxfun; j++)
            {
                if( ((popFit[i][j] < pbestsFit[i][j]) && (optimization == 0)) || ((popFit[i][j] > pbestsFit[i][j]) && (optimization == 1)))
                    sum += 1;
            }
            if (sum == maxfun)
            {
                better = 0;
            }
            //&& pContribution[i] < bpContribution[i])
            else
            {
                if (sum == 0 )//&& pContribution[i] > bpContribution[i])
                    better = 1;
                else
                    better = RandomInt(0, 1);
            }
            if (better == 0)
            {
                for(j = 0; j < maxfun; j++)
                    pbestsFit[i][j] = popFit[i][j];
                for(j = 0; j < maxvar; j++)
                    pbestsVar[i][j] = popVar[i][j];
            }
        }

    /*for(i = 0, n= 0; i < nondomCtr; i++)
        for(j = 0; j < maxfun; j++)
        {
            points[n] = pbestsFit[i][j];
            n++;
        }
    hypeIndicator(bpContribution,popsize,referencePoints,param_k,points,maxfun);*/
}

void save_results(char *archiveName)  /* Escribir resultados en el archivo */
{
	unsigned int i, j;
	FILE *fp;
    /* Abrir archivo para escritura */
	fp = fopen(archiveName, "w");
	for(i = 0; i < nondomCtr; i++)
	{
		for(j = 0; j < maxfun; j++)
			fprintf(fp, "%e ", archiveFit[i][j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

#include "test-fun.h"
/* Fin */
