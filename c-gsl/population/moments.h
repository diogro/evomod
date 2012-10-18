#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include "pop.h"


void population_moments (Population *pop);
void population_phenotype (const gsl_rng * r, Population * pop);
void population_print_moments (const Population *pop, FILE *stream);
void population_write_moments (const Population * pop, FILE * phenotype, FILE * g_corr, FILE * p_corr, FILE * g_var, FILE * p_var, FILE * h_var);
