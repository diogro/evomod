#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include "pop.h"

void population_alloc (const int n_e, const int p, const int m,
    const int burn_in, const int stabilizing, const int selective, const double mu, const double mu_b,
    const double sigma, const double v_e, const gsl_vector *theta,
    const gsl_matrix *omega, Population *pop);

void population_free (Population *pop);

void population_copy (Population *dest, const Population *src);

void population_replicate (Population *dest, const Population *src);

void population_fprintf (const Population *pop, FILE *stream);

void population_fscanf (Population *pop, FILE *stream);

void population_theta_read (Population * pop, FILE *stream);

void population_omega_read (Population * pop, FILE *stream);
