typedef struct {
   double **genes; 
   double **fen;
   double **fen_ad;
   double **Corr;
   double **CorrG;
   double **P;
   double **G;
   double **w;
   double ***B;
   double *fen_med;
   double *ZTarget;
   double sigma; // tamanho das mutacoes
   int Ne; // tamanho da populacao
   int tr; // tracos
   int gp; // genes por traco
   int m; // numero de alelos: diplodia * (gp * tr + gm * (tr*tr-tr)/2 )
   int burn_in;
   int seletiva;
   double mu; //taxa de mutacao
   double mu_B; //taxa de mutacao da matriz B
   struct Population *left;
   struct Population *right;
} Population;
