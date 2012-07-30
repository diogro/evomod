CFLAGS=-O3 -Wall
LDFLAGS=-lpthread -g -debug

MinForModule.exe: main.o funcoes_globais.o 
	    icc -o MinForModule.exe main.o funcoes_globais.o $(LDFLAGS) $(CFLAGS) 

burnin: main-burnin.o funcoes_globais.o 
	    icc -o burnin.exe main-burnin.o funcoes_globais.o $(LDFLAGS) $(CFLAGS) 

main.o: main.c funcoes_globais.o 
	icc -c main.c $(LDFLAGS) $(CFLAGS)

main-burnin.o: main-burnin.c funcoes_globais.o 
	icc -c main-burnin.c $(LDFLAGS) $(CFLAGS)

funcoes_globais.o :
	    icc -c funcoes_globais.c $(LDFLAGS) $(CFLAGS)

clean:
	rm MinForModule.exe funcoes_globais.o main.o burnin.exe main-burnin.o wbar.txt

cleanall :
	rm MinForModule.exe funcoes_globais.o main.o  CoeficienteIntegracao correlacao.txt desviotracos.txt PFinal.txt MatCorrFinal.txt matCorr.txt pcorrelacao pdesviotracos ptracos pvarP pvarG pvarH tracos.txt var_alelos.txt varP.txt varG.txt varH.txt pmutacao mutacao.txt burnin.exe main-burnin.o pcorrelacaoG correlacaoG.txt
