#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "pop.h"

#define FS 1. //forca de selecao
#define classes 10 // Classes nos histogramas
#define BI_Branch_Interval 1428
#define SEL_Branch_Interval 0
#define MAXTR 100
#define EPS 10e-19

double **genes2;
double ***B2;
int individuo, gerac;

double normal(double mu, double sig)
{
   double x1, x2, w, y1;
   
   do {
      x1 = 2.0 * ((double)(rand()) / (RAND_MAX + 1.)) - 1.0;
      x2 = 2.0 * ((double)(rand()) / (RAND_MAX + 1.)) - 1.0;
      w = x1 * x1 + x2 * x2;
   } while ( w >= 1.0 );

   w = sqrt( (-2.0 * log( w ) ) / w );
   y1 = x1 * w;

   return( mu + y1 * sig );
}

void extend (double **mat, double **v, double *d, int n)
{
   int i, j, k;
   double **aux;

   aux = (double **) malloc (n*sizeof(double *));
   for (j = 0; j < n; j++)
      aux[j] = (double *) malloc (n*sizeof(double));
   
   for ( i = 0; i < n; i++)
      if(d[i] < 0.0) d[i] = d[i-1];

   for ( i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         mat[i][j] = 0.0;
      
   for ( i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         aux[i][j] = d[i]*v[j][i]; 

   for ( i = 0; i < n; i++)
      for (j = 0; j < n; j++)
         for (k = 0; k < n; k++)
            mat[i][j] += v[i][k]*aux[k][j];
   free(aux);
}

void eigsrt (double *d, double ***v, int n)
{
   int i, k, j;
   double p;
   for (i = 0; i<n-1; i++){
      p = d[k=i];
      for (j = i; j<n; j++)
         if (d[j] >= p) p=d[k=j];
      if( k!=i ){
         d[k] = d[i];
         d[i] = p;
         for ( j = 0; j<n; j++){
            p = (*v)[j][i];
            (*v)[j][i] =(*v)[j][k];
            (*v)[j][k] = p;
         }
      }
   }
}

void rot (double **a, double s, double tau, int i, int j, int k, int l)
{
   double g = a[i][j];
   double h = a[k][l];
   a[i][j] = g - s*(h+g*tau);
   a[k][l] = h + s*(g-h*tau);
}

void jacobi (double **a, double **v, double *d, int n)
{
   int nrot = 0, i, j, ip, iq;
   double tresh, theta, tau, t, sm, s, h, g, c;
   double *b, *z;

   z = (double *) malloc (n*sizeof(double));
   b = (double *) malloc (n*sizeof(double));

   for (ip = 0; ip < n; ip++){
      for (iq = 0; iq < n; iq++)
         v[ip][iq] = 0.;
   }
   for (ip = 0; ip < n; ip++){
      v[ip][ip] = 1.;
      d[ip] = a[ip][ip];
      b[ip] = d[ip];
      z[ip] = 0.;
   }
   for ( i = 1; i <= 50; i++) {
      sm = 0.0;
      for (ip = 0; ip < n-1; ip++){
         for (iq = ip+1; iq<n; iq++)
            sm += fabs(a[ip][iq]);
      }
      if (sm < EPS){
         eigsrt (d, &v, n);
         free(b);
         free(z);
         return;
      }
      if (i < 4)
         tresh = 0.2*sm/(n*n);
      else
         tresh = 0.0;
      for (ip = 0; ip < n-1; ip++){
         for (iq = ip+1; iq<n; iq++){
            g = 100.0*fabs(a[ip][iq]);
            if ( i > 4 && g <= EPS*fabs(d[ip]) && g <= EPS*fabs(d[iq]))
               a[ip][iq] = 0.0;
            else if (fabs(a[ip][iq]) > tresh){
               h = d[iq]-d[ip];
               if (g <= EPS*fabs(h))
                  t = (a[ip][iq])/h;
               else {
                  theta = 0.5*h/(a[ip][iq]);
                  t = 1.0/(fabs(theta) + sqrt(1.0+theta*theta));
                  if (theta < 0.0) t = -t;
               }
               c = 1.0/sqrt(1+t*t);
               s = t*c;
               tau = s/(1.0+c);
               h = t*a[ip][iq];
               z[ip] -= h;
               z[iq] += h;
               d[ip] -= h;
               d[iq] += h;
               a[ip][iq] = 0.0;
               for ( j = 0; j < ip; j++)
                  rot(a, s, tau, j, ip, j, iq);
               for ( j = ip+1; j < iq; j++)
                  rot(a, s, tau, ip, j, j, iq);
               for ( j = iq+1; j < n; j++)
                  rot(a, s, tau, ip, j, iq, j);
               for ( j = 0; j < n; j++)
                  rot(v, s, tau, j, ip, j, iq);
               ++nrot;
            }
         }
      }
      for (ip=0; ip<n; ip++){
         b[ip] += z[ip];
         d[ip] = b[ip];
         z[ip] = 0.0;
      }
   }
   printf("\nNÃO COPNVERIGUEWD QARARRO\n");
   free(b);
   free(z);
}

void bend (double **mat, int tr)
{
   int j;
   double **v, *d;
   FILE *error;
   
   v = (double **) malloc (tr*sizeof(double *));
   for (j = 0; j < tr; j++)
      v[j] = (double *) malloc (tr*sizeof(double));

   d = (double *) malloc (tr*sizeof(double));

   jacobi (mat, v, d, tr);
   extend(mat, v, d, tr);
   error = fopen("./std.error.txt", "a");
   for (j = 0; j < tr; j++)
      fprintf(error, "%lf ", d[j]);
   fprintf(error, "\n");
   free(v);
   free(d);
   fclose(error);
}

void cholesky (double **orig, double **chol, int tr, int dept)
{
   int i, j, k;

   FILE *error;
start:
   for (i = 0; i < tr; i++) {
      chol[i][i] = orig[i][i];
      for (k = 0; k < i; k++)
         chol[i][i] -= chol[k][i]*chol[k][i];
      if (chol[i][i] <= 0) {
         bend(orig, tr);
         error = fopen("./std.error.txt", "a");
         fprintf(error, "\nERROR: non-positive definite matrix!\n");
         dept++;
         fprintf(error, "\nproblem from %d %f %d\n", i, chol[i][i], dept);
         fclose(error);
         goto start;
      }
      chol[i][i] = sqrt(chol[i][i]);

      for (j = i+1; j < tr; j++) {
         chol[i][j] = orig[i][j];
         for (k = 0; k < i; k++)
            chol[i][j] -= chol[k][i]*chol[k][j];
         chol[i][j] /= chol[i][i];
         chol[j][i] = chol[i][j];
      }
   }
}

void gauss_mutation(double **M, double *mutation, int tr)
{
   int i, j;
   double *b, **chol;

   b = (double *) malloc (tr*sizeof(double));

   chol = (double **) malloc (tr*sizeof(double *));
   for (j = 0; j < tr; ++j) 
      chol[j] = (double *) malloc (tr*sizeof(double));

   cholesky (M, chol, tr, 0);
   for (i = 0; i < tr; i++)
      b[i] = normal(0., 1.);

   for ( i = 0; i < tr; i++)
   {
      mutation[i] = 0.;
      for (j = 0; j <=i; j++)
         mutation[i] += chol[i][j]*b[j];
   }
   free(b);
   free(chol);
}

void mutacao(Population *Pop)
{
   int i, j, k;
   int Ne = Pop->Ne;
   int tr = Pop->tr; 
   int m = Pop->m;
   double rand01;

   for (k = 0; k < Ne; k++)
   {
      for (j = 0; j < m; j++)
      {
         rand01 = (double)(rand())/(RAND_MAX);
         if ( rand01 < Pop->mu ) {Pop->genes[k][j] += normal(0, Pop->sigma);}
         for ( i = 0; i < tr; i++)
         {
            rand01 = (double)(rand())/(RAND_MAX);
            if ( rand01 < Pop->mu_B ) {Pop->B[k][j][i] = (1-Pop->B[k][j][i]);}
         }

      }
   }
}

void FenotipoMedio (Population *Pop)
{
   int k, j;

   for ( j = 0; j < Pop->tr; j++)
   {
      Pop->fen_med[j] = 0.;
      for ( k = 0;  k < Pop->Ne; k++)
         Pop->fen_med[j] += Pop->fen[k][j];
      Pop->fen_med[j] /= Pop->Ne;
   }
}

void fen_pop(Population *Pop)
{
   int i, j, k, l, aux;
   int Ne = Pop->Ne;
   int tr = Pop->tr; 
   int m = Pop->m;
   double exp_aux;

   for (k = 0; k < Ne; k++)
      for (j = 0; j < tr; j++){
         Pop->fen[k][j] = 0;
         Pop->fen_ad[k][j] = 0;
      }

   for (k = 0; k < Ne; k++)
   {

      for ( i = 0; i < tr; i++){
         for ( j = 0; j < m; j++){
            Pop->fen_ad[k][i] += Pop->genes[k][j]*Pop->B[k][j][i];
         }
      }
      for ( j = 0; j < tr; j++)
         Pop->fen[k][j] = Pop->fen_ad[k][j] + normal(0., 0.8);
   }
   FenotipoMedio (Pop);
}

void mat_corr_e_P (Population *Pop)
{
   int i, j, k, l;
   double *x, *y;
   double correlation, covariance;

   x = (double *) malloc (Pop->Ne*sizeof(double));
   y = (double *) malloc (Pop->Ne*sizeof(double));

   double yt, xt;
   double syy = 0.0, sxy = 0.0, sxx = 0.0, ay = 0.0, ax = 0.0;
   
   for (i = 0; i < Pop->tr; i++){
      for (j = 0; j <= i; j++){
         syy = 0.0; 
         sxy = 0.0;
         sxx = 0.0;
         ay = 0.0; 
         ax = 0.0;
         for(k = 0; k < Pop->Ne; k++){
            x[k] = Pop->fen[k][i];
            y[k] = Pop->fen[k][j];
            ax += x[k];
            ay += y[k];
         }
         ax /= Pop->Ne;
         ay /= Pop->Ne;
         for(k = 0; k < Pop->Ne; k++){
            xt = x[k] - ax;
            yt = y[k] - ay;
            sxx += xt*xt;
            syy += yt*yt;
            sxy += xt*yt;
         }
         correlation = sxy/sqrt(sxx*syy); 
         covariance = sxy/Pop->Ne;
         Pop->Corr[i][j] = correlation;
         Pop->Corr[j][i] = correlation;
         Pop->P[i][j] = covariance;   
         Pop->P[j][i] = covariance;
      }
   }
   for (i = 0; i < Pop->tr; i++){
      for (j = 0; j <= i; j++){
         syy = 0.0; 
         sxy = 0.0;
         sxx = 0.0;
         ay = 0.0; 
         ax = 0.0;
         for(k = 0; k < Pop->Ne; k++){
            x[k] = Pop->fen_ad[k][i];
            y[k] = Pop->fen_ad[k][j];
            ax += x[k];
            ay += y[k];
         }
         ax /= Pop->Ne;
         ay /= Pop->Ne;
         for(k = 0; k < Pop->Ne; k++){
            xt = x[k] - ax;
            yt = y[k] - ay;
            sxx += xt*xt;
            syy += yt*yt;
            sxy += xt*yt;
         }
         correlation = sxy/sqrt(sxx*syy); 
         covariance = sxy/Pop->Ne;
         Pop->CorrG[i][j] = correlation;
         Pop->CorrG[j][i] = correlation;
         Pop->G[i][j] = covariance;   
         Pop->G[j][i] = covariance;
      }
   }
   free(x);
   free(y);
}

double solve(double *left, double **chol, double *right, int tr)
{
   int i, k;
   double result = 0.;

   for (i = 0; i < tr; i++) {
      for (k = 0; k < i; k++) {
         right[i] -= right[k]*chol[k][i];
      }
      right[i] /= chol[i][i];
   }
   
   for (i = tr-1; i >= 0; i--){
      for (k = i+1; k < tr; k++) {
         right[i] -= right[k]*chol[i][k];
      }
      right[i] /= chol[i][i];
   }
   
   for (i = 0; i < tr; i++){ 
      result += right[i]*left[i];
   }

   return result/tr;
}

double selecao_correlacao (double *fen_ind, double *ZTarget, double **w_Chol, int Ne, int tr, int gp)
{
   int i, j;
   double *left_DZ, *DZ, aux;

   DZ = (double *) malloc (tr*sizeof(double));
   left_DZ = (double *) malloc (tr*sizeof(double));

   for ( i = 0; i < tr; i++){
      DZ[i] = fen_ind[i] - ZTarget[i];
      left_DZ[i] = FS*(DZ[i]);
   }

   aux = exp (-0.5 * solve(left_DZ, w_Chol, DZ, tr));

   free(left_DZ);
   free(DZ);
   return aux;
}

void rank_pop (Population *Pop, double *ind_w)
{
   /* Calcula o fitness de cada individuo, depois cria um vetor com
    *     *     * fitness cumulativo de cada individuo */

      int i, j, k; 
      double **chol, aux, aux2;
      FILE *W_bar;

      if(gerac == 0)
         W_bar = fopen("wbar.txt", "w");
      else
         W_bar = fopen("wbar.txt", "a");

      chol = (double **) malloc (Pop->tr*sizeof(double *));
      for (j = 0; j < Pop->tr; ++j) 
         chol[j] = (double *) malloc (Pop->tr*sizeof(double));

      cholesky (Pop->w, chol, Pop->tr, 0);

      aux = 0;
      aux2 = 0;
      for (k = 0; k < Pop->Ne; k++){
         ind_w[k] = selecao_correlacao(Pop->fen[k], Pop->ZTarget, chol, Pop->Ne, Pop->tr, Pop->gp);
         aux += ind_w[k];
         aux2 += ind_w[k]*ind_w[k];
      }

      fprintf(W_bar, "%d %lf %lf\n", gerac, aux/Pop->Ne, sqrt(aux2/Pop->Ne - aux/Pop->Ne*aux/Pop->Ne) );

      ind_w[0] /= aux;
      for (k = 1; k < Pop->Ne; k++){
         ind_w[k] /= aux;
         ind_w[k] += ind_w[k-1];
      }
      free(chol);
      fclose(W_bar);
}

void gera_individuos (Population *Pop, double **genes2, double ***B2)
{
   int i, j, k, l;
   int individuo1, individuo2, alelo1, alelo2;
   int Ne = Pop->Ne;
   int gp = Pop->gp;
   int tr = Pop->tr; 
   int m = Pop->m;
   double rand01, rand02, *ind_w;

   if (gerac <= Pop->burn_in){
      for ( k = 0; k < Pop->Ne; k++) { 
         individuo1 = (int)(Pop->Ne*(double)(rand())/(RAND_MAX));
         individuo2 = (int)(Pop->Ne*(double)(rand())/(RAND_MAX));

         for (j = 0; j < m; j+=2)
         {
            alelo1 = j + (int) (2*(double)(rand())/(RAND_MAX));
            alelo2 = j + (int) (2*(double)(rand())/(RAND_MAX));
            genes2[k][j] = Pop->genes[individuo1][alelo1];
            genes2[k][j+1] = Pop->genes[individuo2][alelo2];
            for ( i = 0; i < tr; i++)
            {
               B2[k][j][i] = Pop->B[individuo1][alelo1][i];
               B2[k][j+1][i] = Pop->B[individuo2][alelo2][i];
            }
         }
      }
   }
   else{
      ind_w = (double *) malloc (Pop->Ne*sizeof(double));
      rank_pop (Pop, ind_w);
      for ( k = 0; k < Pop->Ne; k++) { 
         rand01 = (double)(rand())/(RAND_MAX); 
         rand02 = (double)(rand())/(RAND_MAX); 
         individuo1 = 0; 
         while(ind_w[individuo1] < rand01){ individuo1++; }
         individuo2 = 0;
         while(ind_w[individuo2] < rand02){ individuo2++; }
         for (j = 0; j < m; j+=2)
         {
            alelo1 = j + (int) (2*(double)(rand())/(RAND_MAX));
            alelo2 = j + (int) (2*(double)(rand())/(RAND_MAX));
            genes2[k][j] = Pop->genes[individuo1][alelo1];
            genes2[k][j+1] = Pop->genes[individuo2][alelo2];
            for ( i = 0; i < tr; i++)
            {
               B2[k][j][i] = Pop->B[individuo1][alelo1][i];
               B2[k][j+1][i] = Pop->B[individuo2][alelo2][i];
            }
         }
      }
      free(ind_w);
   }
}

void cruzamento(Population *Pop, int geracao)
{
   int i, j, k;

   gerac = geracao;

   gera_individuos (Pop, genes2, B2);
   
   for (k = 0; k < Pop->Ne; k++)
      for (j = 0; j < Pop->m; j++){
         Pop->genes[k][j] = genes2[k][j];
         for ( i = 0; i < Pop->tr; i++)
            Pop->B[k][j][i] = B2[k][j][i];
      }
   
   fen_pop (Pop);
   mat_corr_e_P (Pop);
}

double pearson (double *x, double *y, unsigned long n)
{
   unsigned long j;
   double yt, xt;
   double syy = 0.0, sxy = 0.0, sxx = 0.0, ay = 0.0, ax = 0.0;

   for(j = 0; j < n; j++){
      ax += x[j];
      ay += y[j];
   }
   ax /= n;
   ay /= n;
   for(j = 0; j < n; j++){
      xt = x[j] - ax;
      yt = y[j] - ay;
      sxx += xt*xt;
      syy += yt*yt;
      sxy += xt*yt;
   }
   return sxy/sqrt(sxx*syy); 
}

void print_pop (Population *Pop, char *arquivo)
{
   int i, j, k; 
   FILE *Out;
   Out =  fopen( arquivo, "w" );

   fprintf(Out, "%d\n", Pop->Ne );
   fprintf(Out, "%d\n", Pop->tr );
   fprintf(Out, "%d\n", Pop->gp );
   fprintf(Out, "%lf\n", Pop->mu);
   fprintf(Out, "%lf\n", Pop->mu_B);
   fprintf(Out, "%lf\n", Pop->sigma);

   for (i = 0; i < Pop->Ne; i++)
      for(j = 0; j < Pop->m; j++)
         fprintf(Out, "%lf ", Pop->genes[i][j]);
   for (k = 0; k < Pop->Ne; k++)
      for (i = 0; i < Pop->m; i++) 
         for(j = 0; j < Pop->tr; j++)
            fprintf(Out, "%lf ", Pop->B[k][i][j]);
}

void read_pop (Population *Pop1, char *arquivo)
{

   int i, j, k; 
   FILE *Out;
   Out =  fopen( arquivo, "r" );
   fscanf(Out, "%d", &Pop1->Ne );
   fscanf(Out, "%d", &Pop1->tr );
   fscanf(Out, "%d", &Pop1->gp );
   fscanf(Out, "%lf", &Pop1->mu);
   fscanf(Out, "%lf", &Pop1->mu_B);
   fscanf(Out, "%lf", &Pop1->sigma);
   
   Pop1->m = 2 * (Pop1->gp * Pop1->tr);
   
   Pop1->genes = (double **) malloc (Pop1->Ne*sizeof(double *));
   for (j = 0; j < Pop1->Ne; ++j) 
      Pop1->genes[j] = (double *) malloc (Pop1->m*sizeof(double));

   Pop1->fen_med = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->fen = (double **) malloc (Pop1->Ne*sizeof(double *));
   for (j = 0; j < Pop1->Ne; ++j) 
      Pop1->fen[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->fen_ad = (double **) malloc (Pop1->Ne*sizeof(double *));
   for (j = 0; j < Pop1->Ne; ++j) 
      Pop1->fen_ad[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->Corr = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->Corr[j] = (double *) malloc (Pop1->tr*sizeof(double));

   Pop1->P = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->P[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->CorrG = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->CorrG[j] = (double *) malloc (Pop1->tr*sizeof(double));

   Pop1->G = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->G[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->B = (double ***) malloc (Pop1->Ne*sizeof(double **));
   for (i = 0; i < Pop1->Ne; ++i) {
      Pop1->B[i] = (double **) malloc (Pop1->m*sizeof(double *));
      for (j = 0; j < Pop1->m; ++j) 
         Pop1->B[i][j] = (double *) malloc (Pop1->tr*sizeof(double));
   }
   
   Pop1->ZTarget = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->w = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->w[j] = (double *) malloc (Pop1->tr*sizeof(double));

   for (i = 0; i < Pop1->Ne; i++)
      for(j = 0; j < Pop1->m; j++)
         fscanf(Out, "%lf", &Pop1->genes[i][j]);
   for (k = 0; k < Pop1->Ne; k++)
      for (i = 0; i < Pop1->m; i++) 
         for(j = 0; j < Pop1->tr; j++)
            fscanf(Out, "%lf", &Pop1->B[k][i][j]);
}

Population * replica_pop (Population *Pop2)
{
   int i, j, k;

   Population *Pop1;
   Pop1 = (Population *) malloc (sizeof(Population));

   Pop1->Ne = Pop2->Ne;
   Pop1->tr = Pop2->tr;
   Pop1->gp = Pop2->gp;
   Pop1->m = Pop2->m;
   Pop1->mu = Pop2->mu;
   Pop1->mu_B = Pop2->mu_B;
   Pop1->sigma = Pop2->sigma;
   Pop1->burn_in = Pop2->burn_in;
   Pop1->seletiva = Pop2->seletiva;
   
   Pop1->genes = (double **) malloc (Pop1->Ne*sizeof(double *));
   for (j = 0; j < Pop1->Ne; ++j) 
      Pop1->genes[j] = (double *) malloc (Pop1->m*sizeof(double));
   
   Pop1->fen_med = (double *) malloc (Pop1->tr*sizeof(double));
  
   Pop1->fen = (double **) malloc (Pop1->Ne*sizeof(double *));
   for (j = 0; j < Pop1->Ne; ++j) 
      Pop1->fen[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->Corr = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->Corr[j] = (double *) malloc (Pop1->tr*sizeof(double));

   Pop1->P = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->P[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->fen_ad = (double **) malloc (Pop1->Ne*sizeof(double *));
   for (j = 0; j < Pop1->Ne; ++j) 
      Pop1->fen_ad[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->CorrG = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->CorrG[j] = (double *) malloc (Pop1->tr*sizeof(double));

   Pop1->G = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->G[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->B = (double ***) malloc (Pop1->Ne*sizeof(double **));
   for (i = 0; i < Pop1->Ne; ++i) {
      Pop1->B[i] = (double **) malloc (Pop1->m*sizeof(double *));
      for (j = 0; j < Pop1->m; ++j) 
         Pop1->B[i][j] = (double *) malloc (Pop1->tr*sizeof(double));
   }
   
   Pop1->ZTarget = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->w = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->w[j] = (double *) malloc (Pop1->tr*sizeof(double));


   for ( i = 0; i < Pop2->Ne; i++)
      for ( j = 0; j < Pop2->m; j++)
         Pop1->genes[i][j] = Pop2->genes[i][j];
   
   for ( i = 0; i < Pop2->tr; i++){
      Pop1->fen_med[i] = Pop2->fen_med[i];
      Pop1->ZTarget[i] = Pop2->ZTarget[i];
   }
   for ( i = 0; i < Pop2->Ne; i++)
      for ( j = 0; j < Pop2->tr; j++){
         Pop1->fen[i][j] = Pop2->fen[i][j];
         Pop1->fen_ad[i][j] = Pop2->fen_ad[i][j];
         for ( k = 0; k < Pop2->m; k++)
            Pop1->B[i][k][j] = Pop2->B[i][k][j];
      }

   for ( i = 0; i < Pop2->tr; i++)
      for ( j = 0; j < Pop2->tr; j++){
         Pop1->P[i][j] = Pop2->P[i][j];
         Pop1->Corr[i][j] = Pop2->Corr[i][j];
         Pop1->G[i][j] = Pop2->G[i][j];
         Pop1->CorrG[i][j] = Pop2->CorrG[i][j];
         Pop1->w[i][j] = Pop2->w[i][j];
      }
   return Pop1;
}

void replace_pop (Population *Pop, Population *Pop_init)
{
   int i, j, k;

   Pop->Ne = Pop_init->Ne;
   Pop->tr = Pop_init->tr;
   Pop->gp = Pop_init->gp;
   Pop->m = Pop_init->m;
   Pop->mu = Pop_init->mu;
   Pop->mu_B = Pop_init->mu_B;
   Pop->sigma = Pop_init->sigma;
   Pop->burn_in = Pop_init->burn_in;
   Pop->seletiva = Pop_init->seletiva;

   for ( i = 0; i < Pop_init->Ne; i++)
      for ( j = 0; j < Pop_init->m; j++)
         Pop->genes[i][j] = Pop_init->genes[i][j];
   
   for ( i = 0; i < Pop_init->tr; i++){
      Pop->fen_med[i] = Pop_init->fen_med[i];
      Pop->ZTarget[i] = Pop_init->ZTarget[i];
   }
   for ( i = 0; i < Pop_init->Ne; i++)
      for ( j = 0; j < Pop_init->tr; j++){
         Pop->fen[i][j] = Pop_init->fen[i][j];
         Pop->fen_ad[i][j] = Pop_init->fen_ad[i][j];
         for ( k = 0; k < Pop_init->m; k++)
            Pop->B[i][k][j] = Pop_init->B[i][k][j];
      }

   for ( i = 0; i < Pop_init->tr; i++)
      for ( j = 0; j < Pop_init->tr; j++){
         Pop->P[i][j] = Pop_init->P[i][j];
         Pop->Corr[i][j] = Pop_init->Corr[i][j];
         Pop->G[i][j] = Pop_init->G[i][j];
         Pop->CorrG[i][j] = Pop_init->CorrG[i][j];
         Pop->w[i][j] = Pop_init->w[i][j];
      }
}

void aloca_pop (Population *Pop1)
{
   int i, j, k;
   
   Pop1->m = 2 * (Pop1->gp * Pop1->tr);
   
   Pop1->genes = (double **) malloc (Pop1->Ne*sizeof(double *));
   for (j = 0; j < Pop1->Ne; ++j)
      Pop1->genes[j] = (double *) malloc (Pop1->m*sizeof(double));
   
   Pop1->fen_med = (double *) malloc (Pop1->tr*sizeof(double));
  
   Pop1->fen = (double **) malloc (Pop1->Ne*sizeof(double *));
   for (j = 0; j < Pop1->Ne; ++j) 
      Pop1->fen[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->Corr = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->Corr[j] = (double *) malloc (Pop1->tr*sizeof(double));

   Pop1->P = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->P[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->fen_ad = (double **) malloc (Pop1->Ne*sizeof(double *));
   for (j = 0; j < Pop1->Ne; ++j) 
      Pop1->fen_ad[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->CorrG = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->CorrG[j] = (double *) malloc (Pop1->tr*sizeof(double));

   Pop1->G = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->G[j] = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->B = (double ***) malloc (Pop1->Ne*sizeof(double **));
   for (i = 0; i < Pop1->Ne; ++i) {
      Pop1->B[i] = (double **) malloc (Pop1->m*sizeof(double *));
      for (j = 0; j < Pop1->m; ++j) 
         Pop1->B[i][j] = (double *) malloc (Pop1->tr*sizeof(double));
   }

   Pop1->ZTarget = (double *) malloc (Pop1->tr*sizeof(double));
   
   Pop1->w = (double **) malloc (Pop1->tr*sizeof(double *));
   for (j = 0; j < Pop1->tr; ++j) 
      Pop1->w[j] = (double *) malloc (Pop1->tr*sizeof(double));


   for ( k = 0; k < Pop1->Ne; k++)
      for ( j = 0; j < Pop1->m; j++)
         Pop1->genes[k][j] = 0;
   
   for ( i = 0; i < Pop1->tr; i++){
      Pop1->fen_med[i] = 0;
      Pop1->ZTarget[i] = 0;
   }

   for ( k = 0; k < Pop1->Ne; k++)
      for ( j = 0; j < Pop1->tr; j++){
         Pop1->fen[k][j] = 0;
         Pop1->fen_ad[k][j] = 0;
      }
   for ( i = 0; i < Pop1->tr; i++)
      for ( j = 0; j < Pop1->tr; j++){
         Pop1->P[i][j] = 0;
         Pop1->Corr[i][j] = 0;
         Pop1->G[i][j] = 0;
         Pop1->CorrG[i][j] = 0;
         Pop1->w[i][j] = 0;
      }
   for ( k = 0; k < Pop1->Ne; k++)
      for ( i = 0; i < Pop1->m; i++)
         for ( j = 0; j < Pop1->tr; j++)
            Pop1->B[k][i][j] = 1;

}

void print_fen (Population *Pop, char *arquivo)
{
   int i, j; 
   FILE *Out;
   Out =  fopen( arquivo, "w" );

   for (i = 0; i < Pop->Ne; i++){
      for(j = 0; j < Pop->tr; j++)
         fprintf(Out, "%lf ", Pop->fen[i][j]);
      fprintf(Out, "\n");
   }
}

void print_fen_ad (Population *Pop, char *arquivo)
{
   int i, j; 
   FILE *Out;
   Out =  fopen( arquivo, "w" );

   for (i = 0; i < Pop->Ne; i++){
      for(j = 0; j < Pop->tr; j++)
         fprintf(Out, "%lf ", Pop->fen_ad[i][j]);
      fprintf(Out, "\n");
   }
}

void gera_arvore(Population *Pop, int n)
{
   int l;
   if ( n == 0)
      return ;
   else { 
      fen_pop (Pop);
      mat_corr_e_P (Pop);
      for (l = 0; l < Pop->burn_in+Pop->seletiva; l++){
         printf("%d\n", l);
         mutacao(Pop);
         cruzamento(Pop, l);
      }
      Pop->left = replica_pop (Pop);
      Pop->right = replica_pop (Pop);
      gera_arvore(Pop->right, n-1);
      gera_arvore(Pop->left, n-1);
   }
}

void print_populations (Population *Pop, int n, int i, int RecDep)
{
   char out_pop[20];
   if(Pop==NULL)
      return ;
   else{
      print_populations(Pop->left, n-1, 0, RecDep);
      scanf("%s", out_pop);
      printf("%s %d\n", out_pop, RecDep-n);
      print_pop(Pop, out_pop);
      scanf("%s", out_pop);
      print_fen(Pop, out_pop);
      print_populations(Pop->right, n-1, 1, RecDep);
   }
}
