/* NSGA-II routine (implementation of the 'main' function) */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
#include <time.h>
# include "global.h"
# include "rand.h"

int nreal;
int nbin;
int nobj;
int ncon;
int popsize;
int function;
double pcross_real;
double pcross_bin;
double pmut_real;
double pmut_bin;
double eta_c;
double eta_m;
int ngen;
int nbinmut;
int nrealmut;
int nbincross;
int nrealcross;
int *nbits;
double *min_realvar;
double *max_realvar;
double *min_binvar;
double *max_binvar;
int bitlength;
clock_t  startTime, endTime;
double duration, clocktime;

int main (int argc, char **argv)
{
    int i,j;
    FILE *fp;
    char name[20];

    population *parent_pop;
    population *child_pop;
    population *mixed_pop;

    strcpy(name,argv[1]);
    function = atoi(argv[2]);

    srand((unsigned)time(NULL));
    seed = (double)rand()/(double)RAND_MAX / 100.0;
   /*seed = 0.02;*/
	
    popsize = 100;
    ngen = atoi(argv[3]);
    printf("%s %d %d \n",name,function,ngen);
    switch(function)
    {
         case 600: /* zdt1 */
            nobj = 2;
			ncon = 0;
			nreal = 30;
			min_realvar = (double *)malloc(nreal*sizeof(double));
            max_realvar = (double *)malloc(nreal*sizeof(double));
            for (i=0; i<nreal; i++)
            {
                min_realvar[i] = 0.0;
                max_realvar[i] = 1.0;
            }
			break;
        case 605: /* zdt2 */
			nobj = 2;
			ncon = 0;
			nreal = 30;
			min_realvar = (double *)malloc(nreal*sizeof(double));
            max_realvar = (double *)malloc(nreal*sizeof(double));
            for (i=0; i<nreal; i++)
            {
                min_realvar[i] = 0.0;
                max_realvar[i] = 1.0;
            }
			break;
        case 610: /* zdt3 */
			nobj = 2;
			ncon = 0;
			nreal = 30;
			min_realvar = (double *)malloc(nreal*sizeof(double));
            max_realvar = (double *)malloc(nreal*sizeof(double));
            for (i=0; i<nreal; i++)
            {
                min_realvar[i] = 0.0;
                max_realvar[i] = 1.0;
            }
			break;
        case 615: /* zdt4 */
			nobj = 2;
			ncon = 0;
			nreal = 10;
			min_realvar = (double *)malloc(nreal*sizeof(double));
            max_realvar = (double *)malloc(nreal*sizeof(double));
            min_realvar[0] = 0.0;
            max_realvar[0] = 1.0;
            for (i=1; i<nreal; i++)
            {
                min_realvar[i] = -5.0;
                max_realvar[i] = 5.0;
            }
			break;
        case 620: /* zdt6 */
			nobj = 2;
			ncon = 0;
			nreal = 30;
			min_realvar = (double *)malloc(nreal*sizeof(double));
            max_realvar = (double *)malloc(nreal*sizeof(double));
            for (i=0; i<nreal; i++)
            {
                min_realvar[i] = 0.0;
                max_realvar[i] = 1.0;
            }
			break;
        case 700: case 705: case 710: case 715: case 720: case 725:/* dtlz1 */
			nobj = 3;
			ncon = 0;
			nreal = 12;
			min_realvar = (double *)malloc(nreal*sizeof(double));
            max_realvar = (double *)malloc(nreal*sizeof(double));
            for (i=0; i<nreal; i++)
            {
                min_realvar[i] = 0.0;
                max_realvar[i] = 1.0;
            }
			break;
        case 730:/* dtlz1 */
			nobj = 3;
			ncon = 0;
			nreal = 22;
			min_realvar = (double *)malloc(nreal*sizeof(double));
            max_realvar = (double *)malloc(nreal*sizeof(double));
            for (i=0; i<nreal; i++)
            {
                min_realvar[i] = 0.0;
                max_realvar[i] = 1.0;
            }
			break;
    }

    pcross_real = 0.9;
    pmut_real = 0.0333;
    eta_c = 15;
    eta_m = 20;
    nbin = 0;

    nbinmut = 0;
    nrealmut = 0;
    nbincross = 0;
    nrealcross = 0;
    parent_pop = (population *)malloc(sizeof(population));
    child_pop = (population *)malloc(sizeof(population));
    mixed_pop = (population *)malloc(sizeof(population));
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (child_pop, popsize);
    allocate_memory_pop (mixed_pop, 2*popsize);
    randomize();
    startTime = clock();
    initialize_pop (parent_pop);
    decode_pop(parent_pop);
    evaluate_pop (parent_pop);
    assign_rank_and_crowding_distance (parent_pop);

    for (i=2; i<=ngen; i++)
    {
        selection (parent_pop, child_pop);
        mutation_pop (child_pop);
        decode_pop(child_pop);
        evaluate_pop(child_pop);
        merge (parent_pop, child_pop, mixed_pop);
        fill_nondominated_sort (mixed_pop, parent_pop);
    }
    fp = fopen(strcat(name,"salida.out"), "w");
	for(i = 0; i < popsize; i++)
	{
		for(j = 0; j < nobj; j++)
		{
			fprintf(fp, "%.6f ", parent_pop->ind[i].obj[j]);
		}
		fprintf(fp, "\n");
	}
	endTime = clock();
	duration = ( endTime - startTime ) / (double)CLOCKS_PER_SEC;
	fprintf(stdout, "%f sec\n", duration);
	fclose(fp);

    if (nreal!=0)
    {
        free (min_realvar);
        free (max_realvar);
    }
    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (child_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2*popsize);
    free (parent_pop);
    free (child_pop);
    free (mixed_pop);
    return (0);
}
