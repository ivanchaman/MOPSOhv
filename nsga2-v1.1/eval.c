/* Routine for evaluating population members  */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Routine to evaluate objective function values and constraints for a population */
void evaluate_pop (population *pop)
{
    int i;
    for (i=0; i<popsize; i++)
    {
        evaluate_ind (&(pop->ind[i]));
    }
    return;
}

/* Routine to evaluate objective function values and constraints for an individual */
void evaluate_ind (individual *ind)
{
    int j;


    /*test_problem (ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);*/
    switch(function)
    {
        case 600: /* zdt1 */
            test_zdt1(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
        case 605: /* zdt2 */
            test_zdt2(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
        case 610: /* zdt3 */
            test_zdt3(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
        case 615: /* zdt4 */
            test_zdt4(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
        case 620: /* zdt6 */
            test_zdt6(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
        case 700: /* zdt6 */
            test_dtlz1(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
        case 705: /* zdt6 */
            test_dtlz2(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
        case 710: /* zdt6 */
            test_dtlz3(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
        case 715: /* zdt6 */
            test_dtlz4(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
        case 720: /* zdt6 */
            test_dtlz5(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
        case 725: /* zdt6 */
            test_dtlz6(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
        case 730: /* zdt6 */
            test_dtlz7(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
			break;
    }

    if (ncon==0)
    {
        ind->constr_violation = 0.0;
    }
    else
    {
        ind->constr_violation = 0.0;
        for (j=0; j<ncon; j++)
        {
            if (ind->constr[j]<0.0)
            {
                ind->constr_violation += ind->constr[j];
            }
        }
    }
    return;
}
