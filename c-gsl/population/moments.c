#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "pop.h"

void population_moments (Population *pop)
{
    int k;

    gsl_vector_set_zero (pop->mean_y);
    gsl_matrix_set_zero (pop->mean_b);
    gsl_vector_set_zero (pop->mean_x);
    gsl_vector_set_zero (pop->mean_z);
    for (k = 0; k < pop->n_e; k++){
        gsl_vector_add (pop->mean_y, pop->y[i]);
        gsl_matrix_add (pop->mean_b, pop->b[i]);
        gsl_vector_add (pop->mean_x, pop->x[i]);
        gsl_vector_add (pop->mean_z, pop->y[i]);
    }
    gsl_vector_scale (pop->mean_y, 1./pop->n_e);
    gsl_matrix_scale (pop->mean_b, 1./pop->n_e);
    gsl_vector_scale (pop->mean_x, 1./pop->n_e);
    gsl_vector_scale (pop->mean_z, 1./pop->n_e);
    covariance_calc (pop->x, pop->n_e, pop->g_matrix, pop->corr_g){
    covariance_calc (pop->z, pop->n_e, pop->p_matrix, pop->corr_p){
}           

void covariance_calc (const gsl_vector ** data, const int n_e, gsl_matrix * cov, gsl_matrix * corr){
    
    int k, i, j;
    double data_i[n_e], data_j[n_e];
    double covariance, correlation;

    for( i = 0; i < pop->p; i++){
        for ( j = 0; j <= i; j++){
            for( k = 0; k < g; k++){
                data_i[k] = gsl_vector_get(data[k], i);
                data_j[k] = gsl_vector_get(data[k], j);
            }
            covariance = gsl_stats_covariance(data_i, sizeof(double), data_j, sizeof(double), n_e);
            correlation = gsl_stats_correlation(data_i, sizeof(double), data_j, sizeof(double), n_e);
            gsl_matrix_set(cov, i, j, covariance);
            gsl_matrix_set(cov, j, i, covariance);
            gsl_matrix_set(corr, i, j, correlation);
            gsl_matrix_set(corr, j, i, correlation);
        }
    }
}
