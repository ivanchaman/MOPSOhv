#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

int dim;    /* number of objectives */
//double bound; /*reference*/

int weaklyDominates( double *point1, double *point2, int no_objectives );
double drand( double from, double to );
void rearrangeIndicesByColumn(double *mat, int rows, int columns, int col,int *ind );
void hypeExactRecursive( double* input_p, int pnts, int dim, int nrOfPnts,int actDim, double* bounds, int* input_pvec, double* fitness,	double* rho, int param_k );
void hypeExact( double* val, int popsize, double *boundsVec,	int param_k, double* points, double* rho);
void hypeSampling( double* val, int popsize, double lowerbound,double upperbound, int nrOfSamples, int param_k, double* points,	double* rho );
void hypeIndicator( double* val, int popsize, double *boundsVec, int param_k, double* points );

double drand( double from, double to )
{
	double j;
	j = from + (double)( (to-from)*rand() / ( RAND_MAX + 1.0) );
	return (j);
}
int weaklyDominates( double *point1, double *point2, int no_objectives )
{
	int better;
	int i = 0;
	better = 1;
	while( i < no_objectives && better )
	{
		better = point1[i] <= point2[i];
		i++;
	}
	return better;
}

void rearrangeIndicesByColumn(double *mat, int rows, int columns, int col,int *ind )
/**
 * Internal function used by hypeExact
 */
{
#define  MAX_LEVELS  300
	int  beg[MAX_LEVELS], end[MAX_LEVELS], i = 0, L, R, swap;
	double pref, pind;
	double ref[rows];
	for( i = 0; i < rows; i++ ) {
		ref[i] = mat[ col + ind[i]*columns ];
	}
	i = 0;

	beg[0] = 0; end[0] = rows;
	while ( i >= 0 ) {
		L = beg[i]; R = end[i]-1;
		if( L < R ) {
			pref = ref[ L ];
			pind = ind[ L ];
			while( L < R ) {
				while( ref[ R ] >= pref && L < R )
					R--;
				if( L < R ) {
					ref[ L ] = ref[ R ];
					ind[ L++] = ind[R];
				}
				while( ref[L] <= pref && L < R )
					L++;
				if( L < R) {
					ref[ R ] = ref[ L ];
					ind[ R--] = ind[L];
				}
			}
			ref[ L ] = pref; ind[L] = pind;
			beg[i+1] = L+1; end[i+1] = end[i];
			end[i++] = L;
			if( end[i] - beg[i] > end[i-1] - beg[i-1] ) {
				swap = beg[i]; beg[i] = beg[i-1]; beg[i-1] = swap;
				swap = end[i]; end[i] = end[i-1]; end[i-1] = swap;
			}
		}
		else {
			i--;
		}
	}
}

void hypeExactRecursive( double* input_p, int pnts, int dim, int nrOfPnts,int actDim, double* bounds, int* input_pvec, double* fitness,	double* rho, int param_k )
/**
 * Internal function used by hypeExact
 */
{
	int i, j;
	double extrusion;
	int pvec[pnts];
	double p[pnts*dim];
	for( i = 0; i < pnts; i++ ) {
		fitness[i] = 0;
		pvec[i] = input_pvec[i];
	}
	for( i = 0; i < pnts*dim; i++ )
		p[i] = input_p[i];

	rearrangeIndicesByColumn( p, nrOfPnts, dim, actDim, pvec );

	for( i = 0; i < nrOfPnts; i++ )
	{
		if( i < nrOfPnts - 1 )
			extrusion = p[ (pvec[i+1])*dim + actDim ] - p[ pvec[i]*dim + actDim ];
		else
			extrusion = bounds[actDim] - p[ pvec[i]*dim + actDim ];

		if( actDim == 0 ) {
			if( i+1 <= param_k )
				for( j = 0; j <= i; j++ ) {
					fitness[ pvec[j] ] = fitness[ pvec[j] ]
					                              + extrusion*rho[ i+1 ];
				}
		}
		else if( extrusion > 0 ) {
			double tmpfit[ pnts ];
			hypeExactRecursive( p, pnts, dim, i+1, actDim-1, bounds, pvec,
					tmpfit, rho, param_k );
			for( j = 0; j < pnts; j++ )
				fitness[j] += extrusion*tmpfit[j];
		}
	}
}

void hypeExact( double* val, int popsize,double *boundsVec, int param_k, double* points, double* rho)
/**
 * Calculating the hypeIndicator
 * \f[ \sum_{i=1}^k \left( \prod_{j=1}^{i-1} \frac{k-j}{|P|-j} \right) \frac{ Leb( H_i(a) ) }{ i } \f]
 */
{
	int i;

	int indices[ popsize ];
	//for( i = 0; i < dim; i++ )
		//boundsVec[i] = bound;
	for( i = 0; i < popsize; i++  )
		indices[i] = i;
    //printf("Calculo Exacto\n");
	/** Recursively calculate the indicator values */
	hypeExactRecursive( points, popsize, dim, popsize, dim-1, boundsVec,indices, val, rho, param_k );
}

void hypeSampling( double* val, int popsize, double lowerbound,double upperbound, int nrOfSamples, int param_k, double* points,	double* rho )
/**
 * Sampling the hypeIndicator
 * \f[ \sum_{i=1}^k \left( \prod_{j=1}^{i-1} \frac{k-j}{|P|-j} \right) \frac{ Leb( H_i(a) ) }{ i } \f]
 *
 * @param[out] val vector of all indicators
 * @param[in] popsize size of the population \f$ |P| \f$
 * @param[in] lowerbound scalar denoting the lower vertex of the sampling box
 * @param[in] upperbound scalar denoting the upper vertex of the sampling box
 * @param[in] nrOfSamples the total number of samples
 * @param[in] param_k the variable \f$ k \f$
 * @param[in] points matrix of all objective values dim*popsize entries
 * @param[in] rho weight coefficients
 * @pre popsize >= 0 && lowerbound <= upperbound && param_k >= 1 &&
 * 		param_k <= popsize
 */
{
	assert(popsize >= 0 );
	assert( lowerbound <= upperbound );
	assert( param_k >= 1 );
	assert( param_k <= popsize );

	int i, s, k;
	int hitstat[ popsize ];
	int domCount;

	double sample[ dim ];
	for( s = 0; s < nrOfSamples; s++ )
	{
		for( k = 0; k < dim; k++ )
			sample[ k ] = drand( lowerbound, upperbound );

		domCount = 0;
		for( i = 0; i < popsize; i++ )
		{
			if( weaklyDominates( points + (i*dim), sample, dim) )
			{
				domCount++;
				if( domCount > param_k )
					break;
				hitstat[i] = 1;
			}
			else
				hitstat[i] = 0;
		}
		if( domCount > 0 && domCount <= param_k )
		{
			for( i = 0; i < popsize; i++ )
				if( hitstat[i] == 1 )
					val[i] += rho[domCount];
		}
	}
	for( i = 0; i < popsize; i++ )
	{
		val[i] = val[i] * pow( (upperbound-lowerbound), dim ) / (double)nrOfSamples;
	}
}

void hypeIndicator( double* val, int popsize, double *boundsVec, int param_k, double* points, int dimension)
/**
 * Determine the hypeIndicator
 * \f[ \sum_{i=1}^k \left( \prod_{j=1}^{i-1} \frac{k-j}{|P|-j} \right) \frac{ Leb( H_i(a) ) }{ i } \f]
 *
 * if nrOfSamples < 0, then do exact calculation, else sample the indicator
 *
 * @param[out] val vector of all indicator values
 * @param[in] popsize size of the population \f$ |P| \f$
 * @param[in] lowerbound scalar denoting the lower vertex of the sampling box
 * @param[in] upperbound scalar denoting the upper vertex of the sampling box
 * @param[in] nrOfSamples the total number of samples or, if negative, flag
 * 		that exact calculation should be used.
 * @param[in] param_k the variable \f$ k \f$
 * @param[in] points matrix of all objective values dim*popsize entries
 * @param[in] rho weight coefficients
 */
{
	int i,j;
	double rho[param_k + 1];
    double lowerbound, upperbound;
    int nrOfSamples;
    dim = dimension;
    //printf("cALCULANDO\n");
	/** Set alpha */
	rho[0] = 0;
	for( i = 1; i <= param_k; i++ )
	{
		rho[i] = 1.0 / (double)i;
		for( j = 1; j <= i-1; j++ )
			rho[i] *= (double)(param_k - j ) / (double)( popsize - j );
      //  printf("%f\n", rho[i]);
	}
	for( i = 0; i < popsize; i++ )
		val[i] = 0.0;

	if( dim < 3 )
	{
	    // printf("\nExacto");
		hypeExact( val, popsize, boundsVec, param_k, points,rho );
	}
	else
	{
	   // printf("\nMuestreo");
        nrOfSamples = 10000;
        //dim = dim -1;
        lowerbound = boundsVec[0];
        upperbound = boundsVec[1];
		hypeSampling( val, popsize, lowerbound, upperbound, nrOfSamples,param_k, points, rho);
	}
}

/*

int main()
{
    double val[3];
    int popsize = 2;
    int param_k = 1;
    double points [6] = {1,1,1,3,3,3};
    double pointsref [3] = {0,6,6};
    hypeIndicator(val,popsize, pointsref, param_k, points,3);
    printf("%f %f %f", val[0],val[1],val[2]);
    return 0;
}

*/
