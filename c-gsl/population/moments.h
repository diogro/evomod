#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include "pop.h"


void population_moments (Population *pop);
void population_phenotype (const gsl_rng * r, Population * pop);
