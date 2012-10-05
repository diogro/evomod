#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "pop.h"

void population_alloc (const int n_e, const int p, const int m,
    const int burn_in, const int selective, const double mu, const double mu_b,
    const double sigma, const double v_e, const gsl_vector *theta,
    const gsl_matrix *omega, Population *pop)
{
    int i;
    pop->n_e = n_e;
    pop->p = p;
    pop->m = 2*m;
    pop->burn_in = burn_in;
    pop->selective = selective;
    pop->mu = mu;
    pop->mu_b = mu_b;
    pop->sigma = sigma;
    pop->v_e = v_e;

    pop->y = (gsl_vector **) malloc (pop->n_e*sizeof(gsl_vector *));
    for (i = 0; i < pop->n_e; i++)
        pop->y[i] = gsl_vector_alloc(pop->m);

    pop->b = (gsl_matrix **) malloc (pop->n_e*sizeof(gsl_matrix *));
    for (i = 0; i < pop->n_e; i++)
        pop->b[i] = gsl_matrix_alloc(pop->p, pop->m);

    pop->x = (gsl_vector **) malloc (pop->n_e*sizeof(gsl_vector *));
    for (i = 0; i < pop->n_e; i++)
        pop->x[i] = gsl_vector_alloc(pop->p);

    pop->z = (gsl_vector **) malloc (pop->n_e*sizeof(gsl_vector *));
    for (i = 0; i < pop->n_e; i++)
        pop->z[i] = gsl_vector_alloc(pop->p);

    pop->mean_y = gsl_vector_alloc(pop->m);
    pop->mean_b = gsl_matrix_alloc(pop->p, pop->m);
    pop->mean_x = gsl_vector_alloc(pop->p);
    pop->mean_z = gsl_vector_alloc(pop->p);
    
    pop->g_matrix = gsl_matrix_alloc(pop->p, pop->p);
    pop->p_matrix = gsl_matrix_alloc(pop->p, pop->p);
    pop->corr_g = gsl_matrix_alloc(pop->p, pop->p);
    pop->corr_p = gsl_matrix_alloc(pop->p, pop->p);

    gsl_vector_memcpy (pop->theta, theta)
    gsl_matrix_memcpy (pop->omega, omega)
}
