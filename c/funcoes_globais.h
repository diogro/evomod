#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pop.h"

#define NUM_THREADS 1
#define classes 10
#define MAXTR 100

pthread_mutex_t mutexPop;

double **genes2;
double ***B2;
int individuo, gerac;

double normal(double mu, double sig);
void extend (double **a, double **v, double *d, int n);
void eigsrt (double *d, double ***v, int n);
void rot (double **a, double s, double tau, int i, int j, int k, int l);
void jacobi (double **a, double **v, double *d, int n);
void bend (double **mat, int tr);
void cholesky (double **orig, double **chol, int tr);
void gauss_mutation(double **M, double *mutation, int tr);
void mutacao( Population *Pop );
void fen_pop( Population *Pop );
void mat_corr_e_P ( Population *Pop );
double solve(double *left, double **chol, double *right, int tr);
double selecao_correlacao (double *genes_ind, double *ZTarget, double **GTarget_Chol, int Ne, int tr, int gp);
void rank_pop (Population *Pop, double *ind_w);
void gera_individuos (Population *Pop, double **genes2, double ***B2);
void cruzamento(Population *Pop, int geracao);
void FenotipoMedio (Population *Pop);
double pearson (double *x, double *y, unsigned long n);
void print_pop (Population *Pop, char *arquivo);
void read_pop (Population *Pop1, char *arquivo);
Population * replica_pop (Population *Pop2);
void aloca_pop (Population *Pop2);
void print_fen (Population *Pop, char *arquivo);
void gera_arvore( Population *Pop, int n);
void print_populations (Population *Pop, int n, int i, int RecDep);
