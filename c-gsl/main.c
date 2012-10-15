#include "population/evolution.h"

int main(){

    int n_e, traits, m, burn_in, selective;
    double mu, mu_b, sigma, v_e;

    FILE * phenotype;
    FILE * p_corr;
    FILE * g_corr;
    FILE * g_var;
    FILE * p_var;
    FILE * h_var;
    FILE * out_population;
    FILE * summary;

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

    n_e = 1000;
    traits = 10;
    m = 500;
    mu = 0.0005;
    mu_b = 0.0001;
    sigma = 0.02;
    v_e = 0.8;

    burn_in = 0;
    selective = 1000;

    gsl_vector * theta = gsl_vector_alloc (traits);
    gsl_matrix * omega = gsl_matrix_alloc (traits, traits);

    population_alloc (n_e, traits, m, burn_in, selective, mu, mu_b, sigma, v_e, theta, omega, pop);

    return 0;
}
