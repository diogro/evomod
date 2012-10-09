#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
  
typedef struct {
    int n_e; // population size
    int p; // number of traits
    int m; // number of (diploid) loci
    int burn_in;  // number of drift generations
    int selective; // number of selective generations
    int current_gen; // current generation
    double mu; // genotipic mutation rate
    double mu_b; // ontogenetic mutatation rate
    double sigma; // mutation variance
    double v_e; // envairomental variance
    double *fitness;
    gsl_vector **y; // genome of population
    gsl_matrix **b; // ontogenetic matrix of population
    gsl_vector **x; // additive effects
    gsl_vector **z; // phenotipic values
    gsl_vector *mean_y; // mean phenotipic values
    gsl_matrix *mean_b; // mean additive effects
    gsl_vector *mean_x; // mean additive effects
    gsl_vector *mean_z; // mean phenotipic values
    gsl_matrix *g_matrix; // additive covariance matrix
    gsl_matrix *p_matrix; // phenotipic covarance matrix
    gsl_matrix *corr_g; // additive correlation matrix
    gsl_matrix *corr_p; // phenotipic correlation matrix
    gsl_vector *theta; // optimal position
    gsl_matrix *omega; // selective surface covariance structure
} Population;
