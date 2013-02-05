#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include "population/evolution.h"

int main(){

    int n_e, traits, m, burn_in, selective, generation, i, j;
    double mu, mu_b, sigma, v_e;
    double *delta_theta;


    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_env_setup();

    /*gsl_set_error_handler_off();*/

    FILE * phenotype;
    FILE * p_corr;
    FILE * g_corr;
    FILE * g_var;
    FILE * p_var;
    FILE * h_var;
    FILE * out_population;
    FILE * summary;

    FILE * omega_file;
    FILE * theta_file;

    Population *pop;
    pop = (Population *) malloc (sizeof(Population));

    phenotype      = fopen("./output/phenotype.dat", "w");
    p_corr         = fopen("./output/p.corr.dat", "w");
    g_corr         = fopen("./output/g.corr.dat", "w");
    g_var          = fopen("./output/g.var.dat", "w");
    p_var          = fopen("./output/p.var.dat", "w");
    h_var          = fopen("./output/h.var.dat", "w");
    out_population = fopen("./output/pop.pop", "w");
    summary        = fopen("./output/pop.summary.dat", "w");
    omega_file     = fopen("./omega.csv", "r");
    theta_file     = fopen("./theta.csv", "r");

    n_e       = 5000;
    traits    = 10;
    m         = 500;
    mu        = 0.0005;
    mu_b      = 0.0001;
    sigma     = 0.02;
    v_e       = 0.8;

    burn_in   = 0;
    selective = 11000;

    delta_theta = (double *) malloc (traits*sizeof(double));

    delta_theta[0] = 0.01;
    delta_theta[1] = 0.01;
    delta_theta[2] = 0.01;
    delta_theta[3] = 0.01;
    delta_theta[4] = 0.01;
    delta_theta[5] = -0.01;
    delta_theta[6] = -0.01;
    delta_theta[7] = -0.01;
    delta_theta[8] = -0.01;
    delta_theta[9] = -0.01;

    gsl_vector * theta = gsl_vector_alloc (traits);
    gsl_matrix * omega = gsl_matrix_alloc (traits, traits);

    population_alloc (n_e, traits, m, burn_in, selective, mu, mu_b, sigma, v_e, theta, omega, pop);
    population_random_init (r, pop);
    population_print_moments (pop, summary);

    population_theta_read (pop, theta_file);
    population_omega_read (pop, omega_file);

    /*for (i = 0; i < 10; i++) {*/
        /*printf ("%g\n", gsl_vector_get(pop->theta, i));*/
    /*}*/
    /*for (i = 0; i < 10; i++) {*/
        /*for (j = 0; j < 10; j++) {*/
            /*printf ("%g\n", gsl_matrix_get(pop->omega, i, j));*/
        /*}*/
    /*}*/

    for (generation = 0; generation < pop->burn_in + pop->selective; generation++){
        printf("%d\n", generation);
        population_next_generation(r, pop);
        if (generation > 1000)
            population_theta_change(pop, delta_theta);
        population_write_moments (pop, phenotype, g_corr, p_corr, g_var, p_var, h_var);
    }
    population_print_moments (pop, summary);

    return 0;
}
