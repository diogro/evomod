CFLAGS=-O3 -Wall -I/usr/local/include
LDFLAGS=-lgsl -lgslblas -static

evomod.exe: main.o evo_pop_functions.o
	icc -o MinForModule.exe main.o evo_pop_functions.o $(LDFLAGS) $(CFLAGS) 

burnin: main-burnin.o evo_pop_functions.o
	icc -o burnin.exe main-burnin.o evo_pop_functions.o $(LDFLAGS) $(CFLAGS) 

main.o: main.c evo_pop_functions.o
	icc -c main.c $(LDFLAGS) $(CFLAGS)

burnin.o: burnin.c evo_pop_functions.o 
	icc -c main-burnin.c $(LDFLAGS) $(CFLAGS)

evo_pop_functions.o : moments_pop_functions.o io_pop_functions.o
	icc -c evo_pop_functions.c $(LDFLAGS) $(CFLAGS)

io_pop_functions.o :
	icc -c io_pop_functions.c $(LDFLAGS) $(CFLAGS)

moments_pop_functions.o :
	icc -c moments_pop_functions.c $(LDFLAGS) $(CFLAGS)

clean:
	rm main.o io_pop_functions.o pop_functions.o

cleanall :
	rm main.o 
