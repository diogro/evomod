#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include "io.h"

void covariance_calc (gsl_vector ** data, const int n_e, const int p, gsl_matrix * cov, gsl_matrix * corr)
{
    int ind, i, j;
    double data_i[n_e], data_j[n_e];
    double covariance, correlation;

    for( i = 0; i < p; i++){
        for ( j = 0; j <= i; j++){
            for( ind = 0; ind < n_e; ind++){
                data_i[ind] = gsl_vector_get(data[ind], i);
                data_j[ind] = gsl_vector_get(data[ind], j);
            }
            covariance = gsl_stats_covariance(data_i, 1, data_j, 1, n_e);
            correlation = gsl_stats_correlation(data_i, 1, data_j, 1, n_e);
            gsl_matrix_set(cov, i, j, covariance);
            gsl_matrix_set(cov, j, i, covariance);
            gsl_matrix_set(corr, i, j, correlation);
            gsl_matrix_set(corr, j, i, correlation);
        }
    }
}

void population_moments (Population *pop)
{
    int ind;
    gsl_vector_set_zero (pop->mean_y);
    gsl_matrix_set_zero (pop->mean_b);
    gsl_vector_set_zero (pop->mean_x);
    gsl_vector_set_zero (pop->mean_z);
    for (ind = 0; ind < pop->n_e; ind++){
        gsl_vector_add (pop->mean_y, pop->y[ind]);
        gsl_matrix_add (pop->mean_b, pop->b[ind]);
        gsl_vector_add (pop->mean_x, pop->x[ind]);
        gsl_vector_add (pop->mean_z, pop->z[ind]);
    }
    gsl_vector_scale (pop->mean_y, 1./pop->n_e);
    gsl_matrix_scale (pop->mean_b, 1./pop->n_e);
    gsl_vector_scale (pop->mean_x, 1./pop->n_e);
    gsl_vector_scale (pop->mean_z, 1./pop->n_e);
    covariance_calc (pop->x, pop->n_e, pop->p, pop->g_matrix, pop->corr_g);
    covariance_calc (pop->z, pop->n_e, pop->p, pop->p_matrix, pop->corr_p);
}

void population_phenotype (const gsl_rng * r, Population * pop)
{
    int ind, trait;
    double new_trait;
    for (ind = 0; ind < pop->n_e; ind++){
        gsl_blas_dgemv(CblasNoTrans, 1.0, pop->b[ind], pop->y[ind], 0.0, pop->x[ind]);
        for (trait = 0; trait < pop->p; trait++){
            new_trait = gsl_vector_get(pop->x[ind], trait) + gsl_ran_gaussian(r, pop->v_e);
            gsl_vector_set(pop->z[ind], trait, new_trait);
        }
    }
    population_moments (pop);
}

void population_print_moments (const Population *pop, FILE *stream)
{
    fprintf(stream, "Geracao %d\ntheta = \n", pop->current_gen);
    gsl_vector_fprintf(stream, pop->theta, "%f");
    fprintf(stream, "\nz = \n");
    gsl_vector_fprintf (stream, pop->mean_z, "%f");
    fprintf(stream, "\nx = \n");
    gsl_vector_fprintf (stream, pop->mean_x, "%f");
    fprintf(stream, "\ncorr P = \n");
    gsl_matrix_fprintf (stream, pop->corr_p, "%f");
    fprintf(stream, "\ncorr G = \n");
    gsl_matrix_fprintf (stream, pop->corr_g, "%f");
    fprintf(stream, "\nP = \n");
    gsl_matrix_fprintf (stream, pop->p_matrix, "%f");
    fprintf(stream, "\nG = \n");
    gsl_matrix_fprintf (stream, pop->g_matrix, "%f");
    fprintf(stream, "\nB = \n");
    gsl_matrix_fprintf (stream, pop->mean_b, "%f");
    fprintf(stream, "\ny = \n");
    gsl_vector_fprintf (stream, pop->mean_y, "%f");
    fprintf(stream, "\n---------------------\n");
}

void matrix_lower_trig (const gsl_matrix * mat, gsl_vector * lower_trig, int p)
{
    int i, j, count = 0;
    gsl_vector_set_zero(lower_trig);
    for (i = 1; i < p; i++) {
        for (j = 0; j < i; j++) {
            gsl_vector_set (lower_trig, count, gsl_matrix_get(mat, i, j));
            count++;
        }
    }
}

void population_write_moments (const Population * pop, FILE * phenotype, FILE * g_corr, FILE * p_corr, FILE * g_var, FILE * p_var, FILE * h_var)
{
    int i, j;
    gsl_vector * lower_trig = gsl_vector_alloc ((pop->p*pop->p-pop->p)/2);
    gsl_vector * heridabilities = gsl_vector_alloc (pop->p);
    gsl_vector_view diag_p, diag_g;

    fprintf(phenotype, "%d\n", pop->current_gen);
    gsl_vector_fprintf (phenotype, pop->mean_z, "%f");
    fprintf(phenotype, "\n");

    fprintf(g_corr, "%d\n", pop->current_gen);
    matrix_lower_trig(pop->corr_g, lower_trig, pop->p);
    gsl_vector_fprintf (g_corr, lower_trig, "%f");
    fprintf(g_corr, "\n");

    fprintf(p_corr, "%d\n", pop->current_gen);
    matrix_lower_trig(pop->corr_p, lower_trig, pop->p);
    gsl_vector_fprintf (p_corr, lower_trig, "%f");
    fprintf(p_corr, "\n");

    fprintf(g_var, "%d\n", pop->current_gen);
    diag_g = gsl_matrix_diagonal (pop->g_matrix);
    gsl_vector_fprintf (g_var, &diag_g.vector, "%f");

    fprintf(p_var, "%d\n", pop->current_gen);
    diag_p = gsl_matrix_diagonal (pop->p_matrix);
    gsl_vector_fprintf (p_var, &diag_p.vector, "%f");

    fprintf(h_var, "%d\n", pop->current_gen);
    gsl_vector_memcpy (heridabilities, &diag_g.vector);
    gsl_vector_div (heridabilities, &diag_p.vector);
    gsl_vector_fprintf (h_var, heridabilities, "%f");

    gsl_vector_free (lower_trig);
    gsl_vector_free (heridabilities);
}
