#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include "moments.h"

void population_random_init (const gsl_rng *r, Population * pop);
void population_next_generation (const gsl_rng *r, Population * pop);
