#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include "population/evolution.h"

#define PATH_LEN 1024

/* Build "<dir>/<file>" into out; abort if it does not fit. */
static void join_path(char *out, const char *dir, const char *file)
{
    if (snprintf(out, PATH_LEN, "%s%s", dir, file) >= PATH_LEN) {
        fprintf(stderr, "Path too long: %s%s\n", dir, file);
        exit(EXIT_FAILURE);
    }
}

int main(){

    int n_e, traits, m, burn_in, stabilizing, selective, generation, i;
    double mu, mu_b, sigma, v_e;
    double *delta_theta;
    int bool_sim_type;
    int bool_print_pop;
    long seed;

    gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
    seed = time (NULL) * getpid();
    gsl_rng_set (r, seed);
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
    FILE * parameters_in;
    FILE * parameters_out;

    FILE * input_pop;

    FILE * omega_file;
    FILE * theta_file;

    Population *pop;
    pop = (Population *) malloc (sizeof(Population));

    char out_folder_name[256];
    printf ("\nOutput folder name\n");
    scanf ("%255s", out_folder_name);
    char out_folder_name_path[PATH_LEN];
    join_path(out_folder_name_path, "./output/", out_folder_name);
    mkdir(out_folder_name_path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    printf ("%s", out_folder_name_path);

    char aux[PATH_LEN];
    join_path(aux, out_folder_name_path, "/phenotype.dat");
    phenotype      = fopen(aux, "w");
    join_path(aux, out_folder_name_path, "/p.corr.dat");
    p_corr         = fopen(aux, "w");
    join_path(aux, out_folder_name_path, "/g.corr.dat");
    g_corr         = fopen(aux, "w");
    join_path(aux, out_folder_name_path, "/g.var.dat");
    g_var          = fopen(aux, "w");
    join_path(aux, out_folder_name_path, "/p.var.dat");
    p_var          = fopen(aux, "w");
    join_path(aux, out_folder_name_path, "/h.var.dat");
    h_var          = fopen(aux, "w");
    join_path(aux, out_folder_name_path, "/pop.summary.dat");
    summary        = fopen(aux, "w");
    join_path(aux, out_folder_name_path, "/pop.parameters.txt");
    parameters_out = fopen(aux, "w");

    omega_file     = fopen("./input/omega.csv", "r");
    theta_file     = fopen("./input/theta.csv", "r");

    char input_parameters[256];
    printf ("\nInput parameters file name\n");
    scanf ("%255s", input_parameters);
    char input_folder[PATH_LEN];
    join_path(input_folder, "./input/", input_parameters);
    printf ("%s", input_parameters);
    /*char input_parameters[50] = "./parameter.input.txt";*/

    parameters_in  = fopen(input_folder, "r");

    fscanf(parameters_in, "%d", &n_e);
    fscanf(parameters_in, "%d", &traits);
    fscanf(parameters_in, "%d", &m);
    fscanf(parameters_in, "%lf", &mu);
    fscanf(parameters_in, "%lf", &mu_b);
    fscanf(parameters_in, "%lf", &sigma);
    fscanf(parameters_in, "%lf", &v_e);
    fscanf(parameters_in, "%d", &burn_in);
    fscanf(parameters_in, "%d", &stabilizing);
    fscanf(parameters_in, "%d", &selective);
    delta_theta = (double *) malloc (traits*sizeof(double));
    for(i = 0; i < traits; i++){
        fscanf(parameters_in, "%lf", &delta_theta[i]);
    }

    fprintf(parameters_out, "N_e = %d\n", n_e);
    fprintf(parameters_out, "n.traits = %d\n", traits);
    fprintf(parameters_out, "n.loci = %d\n", m);
    fprintf(parameters_out, "mu = %lf\n", mu);
    fprintf(parameters_out, "mu_B = %lf\n", mu_b);
    fprintf(parameters_out, "sigma = %lf\n", sigma);
    fprintf(parameters_out, "v_e = %lf\n", v_e);
    fprintf(parameters_out, "burn_out = %d\n", burn_in);
    fprintf(parameters_out, "stabilizing = %d\n", stabilizing);
    fprintf(parameters_out, "selective = %d\nDelta theta = ", selective);
    for(i = 0; i < traits; i++){
        fprintf(parameters_out, "%lf ", delta_theta[i]);
    }
    fprintf(parameters_out, "\nseed = %ld\n", seed);

    gsl_vector * theta = gsl_vector_alloc (traits);
    gsl_matrix * omega = gsl_matrix_alloc (traits, traits);

    population_alloc (n_e, traits, m, burn_in, stabilizing, selective, mu, mu_b, sigma, v_e, theta, omega, pop);

    population_theta_read (pop, theta_file);
    population_omega_read (pop, omega_file);

    printf ("\nType 0 for new population, 1 to read from file\n");
    scanf ("%d", &bool_sim_type);

    printf ("\nPrint final pop?\n");
    scanf ("%d", &bool_print_pop);

    if (bool_print_pop) {
        join_path(aux, out_folder_name_path, "/pop.pop");
        out_population = fopen(aux, "w");
    }

    if (!bool_sim_type){
        population_random_init (r, pop);
    }
    else{
        char input_file_pop[PATH_LEN];
        printf ("\nInput population file name\n");
        scanf ("%1023s" ,input_file_pop);
        printf ("%s", input_file_pop);
        input_pop  = fopen(input_file_pop, "r");
        population_fscanf (pop, input_pop);
    }

    population_moments (pop);
    population_print_moments (pop, summary);
    pop->burn_in += pop->current_gen;
    pop->stabilizing += pop->burn_in;
    pop->selective += pop->stabilizing;

	for (generation = 0; generation < burn_in + stabilizing + selective; generation++){
		printf("%d\n", generation);
		population_next_generation(r, pop);
		if(generation > burn_in + stabilizing){
			population_theta_change(pop, delta_theta);
		}
		population_write_moments (pop, phenotype, g_corr, p_corr, g_var, p_var, h_var);
	}
    population_print_moments (pop, summary);
    if (bool_print_pop)
        population_fprintf (pop, out_population);

    return 0;
}
