#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "funcoes_globais.h"

int main(void){

   int i, j, k, l;
   
   double *CorrAnterior; // Para calculo de correlacao de matriz
   double *CorrAtual; 
   double **Mut;
   double aux;

   int numero_variancias = 0; // Quantidade de razoes entre variancias alelicas e fenotipicas (um monte...)
   
   FILE *Correlacao;
   FILE *CorrelacaoG;
   FILE *VarP;
   FILE *VarG;
   FILE *VarH;
   FILE *Tracos;
   FILE *CorrelacaoMatriz;
   FILE *DesvioTracos;
   FILE *MatCorrFinal;
   FILE *PFinal;
   FILE *Coef_inte;
   FILE *Ontogenia;
   FILE *G_Target;
   FILE *PCorrelacao;
   FILE *PCorrelacaoG;
   FILE *PVarP;
   FILE *PVarG;
   FILE *PVarH;
   FILE *PTracos;
   FILE *PDesvioTracos;
   
   G_Target = fopen( "/home/diogro/MainProject/Modularidade/MatrizB/Direcional/MinForModule/G-target.csv", "r" );

   Correlacao = fopen( "correlacao.txt", "w" );
   CorrelacaoG = fopen( "correlacaoG.txt", "w" );
   VarP = fopen( "varP.txt", "w" );
   VarG = fopen( "varG.txt", "w" );
   VarH = fopen( "varH.txt", "w" );
   Tracos = fopen( "tracos.txt", "w" );
   CorrelacaoMatriz = fopen( "matCorr.txt", "w" );
   DesvioTracos = fopen( "desviotracos.txt", "w" );
   MatCorrFinal = fopen("MatCorrFinal.txt", "w");
   PFinal = fopen("PFinal.txt", "w");
   Coef_inte = fopen("CoeficienteIntegracao", "w");

   PCorrelacao = fopen( "pcorrelacao", "w" );
   PCorrelacaoG = fopen( "pcorrelacaoG", "w" );
   PVarP = fopen( "pvarP", "w" );
   PVarG = fopen( "pvarG", "w" );
   PVarH = fopen( "pvarH", "w" );
   PTracos = fopen( "ptracos", "w" );
   PDesvioTracos = fopen( "pdesviotracos", "w" );

   srand((unsigned int)time(NULL));

   scanf("%lf", &aux);
   printf("Delta Total = %lf\n", aux);

   Population *Pop;
   Pop = (Population *) malloc (sizeof(Population));
    read_pop(Pop, "/home/diogro/MainProject/Modularidade/MatrizB/Direcional/MinForModule/Pop-Burnin.pop");

   Pop->burn_in = 0;
   Pop->seletiva = 500;

   /*Leitura da matriz de w_ij*/
   for (i = 0; i < Pop->tr; i++){
      for (j = 0; j < Pop->tr; j++){
         l = fscanf(G_Target, "%lf", &Pop->w[i][j]);
         printf ("%lf ", Pop->w[i][j]);
      }
      printf("\n");
   }
   
   for (i = 1; i < Pop->tr; i++)
      for (j = 0; j < i; j++){
         Pop->w[i][j] *= sqrt(Pop->w[i][i]*Pop->w[j][j]);
         Pop->w[j][i] = Pop->w[i][j];
      }

   for (i = 0; i < Pop->tr; i++){
      for (j = 0; j < Pop->tr; j++){
         printf ("%lf ", Pop->w[i][j]);
      }
      printf("\n");
   }

   for ( i = 0; i < Pop->tr; i++){
      if (i == 0)
         fprintf(PDesvioTracos, "plot 'desviotracos.txt' using 1:%d title '%d'", i+2, i+1);
      else
         fprintf(PDesvioTracos, ", 'desviotracos.txt' using 1:%d title '%d'", i+2, i+1);
   }
   fprintf(PDesvioTracos, "\n");

   for ( i = 0; i < Pop->tr; i++){
      if (i == 0)
         fprintf(PTracos, "plot 'tracos.txt' using 1:%d title '%d'", i+2, i+1);
      else
         fprintf(PTracos, ", 'tracos.txt' using 1:%d title '%d'", i+2, i+1);
   }
   fprintf(PTracos, "\n");

   for ( i = 0; i < Pop->tr; i++){
      if (i == 0){
         fprintf(PVarP, "plot 'varP.txt' using 1:%d title '%d'", i+2, i+1);
         fprintf(PVarG, "plot 'varG.txt' using 1:%d title '%d'", i+2, i+1);
         fprintf(PVarH, "plot 'varH.txt' using 1:%d title '%d'", i+2, i+1);
      }
      else{
         fprintf(PVarP, ", 'varP.txt' using 1:%d title '%d'", i+2, i+1);
         fprintf(PVarG, ", 'varG.txt' using 1:%d title '%d'", i+2, i+1);
         fprintf(PVarH, ", 'varH.txt' using 1:%d title '%d'", i+2, i+1);
      }
   }
   fprintf(PVarP, "\n");
   fprintf(PVarG, "\n");
   fprintf(PVarH, "\n");

   for ( i = 0; i < (Pop->tr*Pop->tr-Pop->tr)/2; i++){
      if (i == 0){
         fprintf(PCorrelacao, "plot 'correlacao.txt' using 1:%d w l title ''", i+2);
         fprintf(PCorrelacaoG, "plot 'correlacaoG.txt' using 1:%d w l title ''", i+2);
      }
      else{
         fprintf(PCorrelacao, ", 'correlacao.txt' using 1:%d w l title ''", i+2);
         fprintf(PCorrelacaoG, ", 'correlacaoG.txt' using 1:%d w l title ''", i+2);
      }
   }
   fprintf(PCorrelacao, "\n");
   fprintf(PCorrelacaoG, "\n");

   /*Alocacao das matrizes auxiliares*/
   CorrAnterior = (double *) malloc (((Pop->tr*Pop->tr-Pop->tr)/2)*sizeof(double));
   CorrAtual = (double *) malloc (((Pop->tr*Pop->tr-Pop->tr)/2)*sizeof(double));
   
   genes2 = (double **) malloc (Pop->Ne*sizeof(double *));
   for (j = 0; j < Pop->Ne; ++j){ 
      genes2[j] = (double *) malloc (Pop->m*sizeof(double));
   }
   B2 = (double ***) malloc (Pop->Ne*sizeof(double **));
   for (i = 0; i < Pop->Ne; ++i) {
      B2[i] = (double **) malloc (Pop->m*sizeof(double *));
      for (j = 0; j < Pop->m; ++j) 
         B2[i][j] = (double *) malloc (Pop->tr*sizeof(double));
   }
   
   /*Laço principal*/
   for ( i = 0; i < Pop->tr; i++){
      Pop->ZTarget[i] = 0.;
   }

   fen_pop (Pop);
   mat_corr_e_P (Pop);
   for (l = 0; l < Pop->burn_in+Pop->seletiva; l++){

      printf("%d\n", l);
   
      for ( i = 0; i < Pop->tr; i++){
         if(i < 5)
            Pop->ZTarget[i] -= aux/10000.;//Pop->fen_med[i];
         else
            Pop->ZTarget[i] += aux/10000.;//Pop->fen_med[i];
      }

      k = 0;
      for ( i = 1; i < Pop->tr; i++)
         for ( j = 0; j < i; j++){
            CorrAnterior[k] = Pop->Corr[i][j];
            k++;
         }
      /* Comandos Efetivos */ 
      
      mutacao(Pop);
      cruzamento(Pop, l);
      
      /* Fim dos Comandos Efetivos */  

      k = 0;
      for ( i = 1; i < Pop->tr; i++)
         for ( j = 0; j < i; j++){
            CorrAtual[k] = Pop->Corr[i][j];
            k++;
         }

      /*Imprime a correlacao entre as matrizes de correlacao de duas
       * geracoes sub-sequentes*/
      if(l > 1)
      fprintf(CorrelacaoMatriz, "%d %lf\n", l, pearson(CorrAtual, CorrAnterior, (Pop->tr*Pop->tr-Pop->tr)/2));

      if(l >= 0)
      {
      /* Imprime a evolucao das correlacoes  */
      fprintf(Correlacao ,"%d", l);
      fprintf(CorrelacaoG ,"%d", l);
      for (i = 1; i < Pop->tr; i++)
         for (j = 0; j < i; j++){
            fprintf(Correlacao ," %lf", Pop->Corr[i][j]);
            fprintf(CorrelacaoG ," %lf", Pop->CorrG[i][j]);
         }
      fprintf(Correlacao ,"\n");
      fprintf(CorrelacaoG ,"\n");
      
      /* Imprime a evolucao das Variancias  */
      fprintf(VarP ,"%d", l);
      fprintf(VarG ,"%d", l);
      fprintf(VarH ,"%d", l);
      for (i = 0; i < Pop->tr; i++){
         fprintf(VarP ," %lf", Pop->P[i][i]);
         fprintf(VarG ," %lf", Pop->G[i][i]);
         fprintf(VarH ," %lf", Pop->G[i][i]/Pop->P[i][i]);
      }
      fprintf(VarP ,"\n");
      fprintf(VarG ,"\n");
      fprintf(VarH ,"\n");

      /* Imprime a evolucao dos Tracos  */
      fprintf(Tracos ,"%d", l);
      for (i = 0; i < Pop->tr; i++){
         fprintf(Tracos ," %lf", Pop->fen_med[i]);
      }
      fprintf(Tracos ,"\n");
      }
   }
   /*Imprime a matriz de correlacao final*/
   for (i = 0; i < Pop->tr; i++){
      for (j = 0; j < Pop->tr; j++)
         fprintf(MatCorrFinal ,"%lf\t", Pop->CorrG[i][j]);
      fprintf(MatCorrFinal,"\n");
   }
   /*Imprime a matriz P final*/
   for (i = 0; i < Pop->tr; i++){
      for (j = 0; j < Pop->tr; j++)
         fprintf(PFinal ,"%lf\t", Pop->P[i][j]);
      fprintf(PFinal,"\n");
   }

   print_pop(Pop, "Pop.pop");
   print_fen(Pop, "Pop.fen");

   fclose(Coef_inte);
   fclose(Correlacao);
   fclose(CorrelacaoG);
   fclose(MatCorrFinal);
   fclose(PFinal);
   fclose(VarP);
   fclose(VarG);
   fclose(VarH);
   fclose(Tracos);
   fclose(CorrelacaoMatriz);
   fclose(DesvioTracos);

   fclose(PCorrelacao);
   fclose(PCorrelacaoG);
   fclose(PVarP);
   fclose(PVarG);
   fclose(PVarH);
   fclose(PTracos);
   fclose(PDesvioTracos);

   return 0;
}
