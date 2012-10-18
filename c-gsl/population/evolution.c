#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_errno.h>
#include "moments.h"

void mutate_ind(const gsl_rng *r, Population *pop, const int ind)
{
    unsigned int mutation_num, mut_i;
    unsigned int *b_mut_idx, *b_pos, bi;
    unsigned int allele, b_rand, b_x, b_y, new_b;
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
        for (bi = 0; bi < pop->m*pop->p; bi++) {
            b_mut_idx[bi] = bi;
        }
        gsl_ran_choose(r, b_pos, mutation_num, b_mut_idx, pop->m*pop->p, sizeof(unsigned int));
        for(mut_i = 0; mut_i < mutation_num; mut_i++) {
           b_rand = gsl_rng_uniform_int(r, 1);
           b_y = b_pos[mut_i] % pop->m;
           b_x = (unsigned int)((double)(b_pos[mut_i]) / (double)(pop->m));
           new_b = (unsigned int)(gsl_matrix_get(pop->b[ind], b_x, b_y) + b_rand) % 2;
           gsl_matrix_set(pop->b[ind], b_x, b_y, new_b);
        }
        free(b_mut_idx);
        free(b_pos);
    }
}

void population_mutate (const gsl_rng *r, Population * pop)
{
    int k;
    for ( k = 0; k < pop->n_e; k++){
        mutate_ind(r, pop, k);
    }
}

double fitness_ind (Population * pop, const int ind, const gsl_matrix * omega_cholesky)
{
    gsl_vector * theta_dist = gsl_vector_alloc (pop->p);
    gsl_vector * theta_dist_t = gsl_vector_alloc (pop->p);
    gsl_vector_memcpy (theta_dist, pop->z[ind]);
    gsl_vector_sub (theta_dist, pop->theta);
    gsl_linalg_cholesky_solve (omega_cholesky, theta_dist, theta_dist_t);
    gsl_blas_ddot (theta_dist, theta_dist_t, &pop->fitness[ind]);
    pop->fitness[ind] = gsl_sf_exp((-1./2.)*(pop->fitness[ind]));
    if(!gsl_finite (pop->fitness[ind]))
        pop->fitness[ind] = 0.;
    gsl_vector_free(theta_dist);
    gsl_vector_free(theta_dist_t);
    return pop->fitness[ind];
}

void population_fitness (Population * pop)
{
    int ind;
    double total_fitness = 0.;
    gsl_matrix * omega_cholesky = gsl_matrix_alloc (pop->p, pop->p);
    gsl_matrix_memcpy (omega_cholesky, pop->omega);
    gsl_linalg_cholesky_decomp (omega_cholesky);
    for (ind = 0; ind < pop->n_e; ind++){
        total_fitness += fitness_ind (pop, ind, omega_cholesky);
    }
    if (total_fitness < 0.000000001){
        for (ind = 0; ind < pop->n_e; ind++){
            pop->fitness[ind] = 1.;  /*TODO: conferir essa parada...*/
        }
    }
    gsl_matrix_free(omega_cholesky);
}

void choose_mates (const gsl_rng *r, Population * pop, int * mates)
{
    unsigned int *n, k, i, ind = 0;
    n = (unsigned int *) malloc(pop->n_e*sizeof(unsigned int));
    population_fitness(pop);
    gsl_ran_multinomial (r, pop->n_e, pop->n_e, pop->fitness, n);
    for ( k = 0; k < pop->n_e; k++){
        if ( n[k] > 0){
            for ( i = 0; i < n[k]; i++){
                mates[ind] = k;
                ind++;
            }
        }
    }
    gsl_ran_shuffle (r, mates, pop->n_e, sizeof(int));
    free(n);
}

void cross_ind (const gsl_rng * r, const int m, const int p,
        const gsl_vector * ind_1_y, const gsl_matrix * ind_1_b,
        const gsl_vector * ind_2_y, const gsl_matrix * ind_2_b,
        gsl_vector * new_ind_y, gsl_matrix * new_ind_b)
{
    int k, locus;
    unsigned int *aleles;
    gsl_vector * b_column  = gsl_vector_alloc(p);
    aleles = (unsigned int *) malloc(m*sizeof(unsigned int));
    for (k = 0; k < m; k++){
        aleles[k] = gsl_ran_bernoulli (r, 0.5);
    }
    for (locus = 0; locus < m/2; locus++){
        gsl_vector_set (new_ind_y, 2 * locus    , gsl_vector_get (ind_1_y, 2 * locus+aleles[2 * locus  ]));
        gsl_vector_set (new_ind_y, 2 * locus + 1, gsl_vector_get (ind_2_y, 2 * locus+aleles[2 * locus + 1]));

        gsl_matrix_get_col (b_column, ind_1_b, 2 * locus+aleles[2 * locus]);
        gsl_matrix_set_col (new_ind_b, 2 * locus, b_column);

        gsl_matrix_get_col (b_column, ind_2_b, 2 * locus+aleles[2 * locus + 1]);
        gsl_matrix_set_col (new_ind_b, 2 * locus + 1, b_column);
    }
    gsl_vector_free(b_column);
    free(aleles);
}

void population_cross (const gsl_rng *r, Population * pop)
{
    int k;

    /* TODO: Deixar esses caras globais? */
    int *dames, * sires;
    gsl_vector ** new_pop_y;
    gsl_matrix ** new_pop_b;
    new_pop_y = (gsl_vector **) malloc (pop->n_e*sizeof(gsl_vector *));
    new_pop_b = (gsl_matrix **) malloc (pop->n_e*sizeof(gsl_matrix *));
    for (k = 0; k < pop->n_e; k++){
        new_pop_y[k] = gsl_vector_alloc(pop->m);
        new_pop_b[k] = gsl_matrix_alloc(pop->p, pop->m);
    }
    sires = (int *) malloc(pop->n_e*sizeof(int));
    dames = (int *) malloc(pop->n_e*sizeof(int));
    /* Ou complica a paralelização? */

    choose_mates(r, pop, sires);
    choose_mates(r, pop, dames);
    for ( k = 0; k < pop->n_e; k++){
        cross_ind (r, pop->m, pop->p,
                   pop->y[sires[k]],  pop->b[sires[k]],
                   pop->y[dames[k]],  pop->b[dames[k]],
                   new_pop_y[k],      new_pop_b[k]);
    }
    for ( k = 0; k < pop->n_e; k++){
        gsl_vector_memcpy(pop->y[k], new_pop_y[k]);
        gsl_vector_free(new_pop_y[k]);
        gsl_matrix_memcpy(pop->b[k], new_pop_b[k]);
        gsl_matrix_free(new_pop_b[k]);
    }
    free(new_pop_y);
    free(new_pop_b);
    free(sires);
    free(dames);
    population_phenotype (r, pop);
}

void population_random_init (const gsl_rng *r, Population * pop)
{
    int i;
    double mu, mu_b;
    for (i = 0; i < pop->n_e; i++){
        gsl_vector_set_zero(pop->y[i]);
        gsl_matrix_set_zero(pop->b[i]);
    }
    pop->current_gen = 0;
    mu               = pop->mu;
    mu_b             = pop->mu_b;
    pop->mu          = 1.;
    pop->mu_b        = 0.5;
    population_mutate(r, pop);
    pop->mu          = mu;
    pop->mu_b        = mu_b;
    population_phenotype(r, pop);
}

void population_next_generation (const gsl_rng *r, Population * pop)
{
    population_mutate(r, pop);
    population_cross(r, pop);
    pop->current_gen++;
}
