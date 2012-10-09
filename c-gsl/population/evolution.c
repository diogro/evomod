#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h> 
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include "pop.h"

void mutate_ind(const gsl_rng *r, Population *pop, const int ind)
{
    unsigned int mutation_num, mut_i;
    unsigned int *b_mut_idx, *b_pos, bi;
    unsigned int allele, b_rand, b_x, b_y;
    double new_allele;

    mutation_num = gsl_ran_binomial(r, pop->mu, pop->m);

    if (mutation_num > 0) {
        for(mut_i = 0; mut_i < mutation_num; mut_i++) {
           allele = gsl_rng_uniform_int(r, pop->m); 
           new_allele = gsl_vector_get(pop->y[ind], allele) + gsl_ran_gaussian(r, pop->sigma); 
           gsl_vector_set(pop->y[ind], allele, new_allele);
        }
    }

    mutation_num = gsl_ran_binomial(r, pop->mu_b, pop->m*pop->p); 

    if (mutation_num > 0) {
        b_mut_idx = (unsigned int *) malloc(pop->m*pop->p*sizeof(unsigned int));
        b_pos = (unsigned int *) malloc(mutation_num*sizeof(unsigned int));
        
        for (bi = 0; bi <= pop->m*pop->p; ++bi) {
            b_mut_idx[bi] = bi;
        }

        gsl_ran_choose(r, b_pos, mutation_num, b_mut_idx, pop->m*pop->p, sizeof(unsigned int));

        for(mut_i = 0; mut_i < mutation_num; mut_i++) {
           b_rand = gsl_rng_uniform_int(r, 1); 
           b_x = b_pos[mut_i] % pop->m;
           b_y = b_pos[mut_i] / pop->m; 
           gsl_matrix_set(pop->b[ind], b_x, b_y, (unsigned int)(gsl_matrix_get(pop->b[ind], b_x, b_y) + b_rand) % 2);
        }
    }
}

void population_mutate (const gsl_rng *r, Population * pop)
{
    int k;
    for ( k = 0; k < pop->n_e; k++){
        mutate_ind(r, pop, k);
    }
}

void fitness_ind (Population * pop, const int ind, const gsl_matrix * omega_cholesky)
{
    gsl_vector * theta_dist = gsl_vector_alloc (pop->p);
    gsl_vector * aux = gsl_vector_alloc (pop->p);
    gsl_vector_memcpy (theta_dist, pop->z[ind]);
    gsl_vector_sub (theta_dist, pop->theta);
    gsl_linalg_cholesky_solve (omega_cholesky, theta_dist, aux);
    gsl_blas_ddot (theta_dist, aux, &pop->fitness[ind]);
    pop->fitness[ind] = gsl_sf_exp((-1./2.)*(pop->fitness[ind]));
}

void population_fitness (Population * pop)
{
    int k;
    gsl_matrix * omega_cholesky = gsl_matrix_alloc (pop->p, pop->p);
    gsl_matrix_memcpy (omega_cholesky, pop->omega);
    gsl_linalg_cholesky_decomp (omega_cholesky);
    for (k = 0; k < pop->n_e; k++) {
        fitness_ind (pop, k, omega_cholesky);
    }
}
