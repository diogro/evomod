#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
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
    pop->fitness = (double *) malloc (pop->n_e*sizeof(double));
    pop->y = (gsl_vector **) malloc (pop->n_e*sizeof(gsl_vector *));
    pop->b = (gsl_matrix **) malloc (pop->n_e*sizeof(gsl_matrix *));
    pop->x = (gsl_vector **) malloc (pop->n_e*sizeof(gsl_vector *));
    pop->z = (gsl_vector **) malloc (pop->n_e*sizeof(gsl_vector *));
    for (i = 0; i < pop->n_e; i++){
        pop->y[i] = gsl_vector_alloc(pop->m);
        pop->b[i] = gsl_matrix_alloc(pop->p, pop->m);
        pop->x[i] = gsl_vector_alloc(pop->p);
        pop->z[i] = gsl_vector_alloc(pop->p);
    }
    pop->mean_y = gsl_vector_alloc(pop->m);
    pop->mean_b = gsl_matrix_alloc(pop->p, pop->m);
    pop->mean_x = gsl_vector_alloc(pop->p);
    pop->mean_z = gsl_vector_alloc(pop->p);
    pop->g_matrix = gsl_matrix_alloc(pop->p, pop->p);
    pop->p_matrix = gsl_matrix_alloc(pop->p, pop->p);
    pop->corr_g = gsl_matrix_alloc(pop->p, pop->p);
    pop->corr_p = gsl_matrix_alloc(pop->p, pop->p);
    pop->theta = gsl_vector_alloc(pop->p);
    pop->omega = gsl_matrix_alloc(pop->p, pop->p);
    gsl_vector_memcpy (pop->theta, theta);
    gsl_matrix_memcpy (pop->omega, omega);
}

void population_free (Population *pop)
{
    int i;
    for (i = 0; i < pop->n_e; i++){
        gsl_vector_free(pop->y[i]);
        gsl_matrix_free(pop->b[i]);
        gsl_vector_free(pop->x[i]);
        gsl_vector_free(pop->z[i]);
    }
    free(pop->y);
    free(pop->b);
    free(pop->x);
    free(pop->z);
    free(pop->fitness);
    gsl_vector_free(pop->mean_y);
    gsl_matrix_free(pop->mean_b);
    gsl_vector_free(pop->mean_x);
    gsl_vector_free(pop->mean_z);
    gsl_matrix_free(pop->g_matrix);
    gsl_matrix_free(pop->p_matrix);
    gsl_matrix_free(pop->corr_g);
    gsl_matrix_free(pop->corr_p);
    gsl_vector_free(pop->theta);
    gsl_matrix_free(pop->omega);
    free(pop);
}

void population_copy (Population *dest, const Population *src)
{
    int i;
    dest->n_e = src->n_e;
    dest->p = src->p;
    dest->m = src->m;
    dest->burn_in = src->burn_in;
    dest->selective = src->selective;
    dest->mu = src->mu;
    dest->mu_b = src->mu_b;
    dest->sigma = src->sigma;
    dest->v_e = src->v_e;
    for (i = 0; i < dest->n_e; i++){
        dest->fitness[i] = src->fitness[i];
        gsl_vector_memcpy(dest->y[i], src->y[i]);
        gsl_matrix_memcpy(dest->b[i], src->b[i]);
        gsl_vector_memcpy(dest->x[i], src->x[i]);
        gsl_vector_memcpy(dest->z[i], src->z[i]);
    }
    gsl_vector_memcpy(dest->mean_y, src->mean_y);
    gsl_matrix_memcpy(dest->mean_b, src->mean_b);  
    gsl_vector_memcpy(dest->mean_x, src->mean_x);   
    gsl_vector_memcpy(dest->mean_z, src->mean_z);   
    gsl_matrix_memcpy(dest->g_matrix, src->g_matrix); 
    gsl_matrix_memcpy(dest->p_matrix, src->p_matrix); 
    gsl_matrix_memcpy(dest->corr_g, src->corr_g);   
    gsl_matrix_memcpy(dest->corr_p, src->corr_p);   
    gsl_vector_memcpy (dest->theta, src->theta);
    gsl_matrix_memcpy (dest->omega, src->omega);
}

void population_replicate (Population *dest, const Population *src)
{
    population_alloc (src->n_e, src->p, src->m, src->burn_in, src->selective,
            src->mu, src->mu_b, src->sigma, src->v_e, src->theta, src->omega,
            dest);
    population_copy (dest, src);
}

void population_fprintf (const Population *pop, FILE *stream)
{
    int i;
    for (i = 0; i < pop->n_e; i++){
        gsl_vector_fprintf(stream, pop->y[i], "%f");
        gsl_matrix_fprintf(stream, pop->b[i], "%f");
    }
}

void population_fscanf (Population *pop, FILE *stream)
{
    int i;
    for (i = 0; i < pop->n_e; i++){
        gsl_vector_fscanf(stream, pop->y[i]);
        gsl_matrix_fscanf(stream, pop->b[i]);
    }
} 
