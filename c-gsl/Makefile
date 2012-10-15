CFLAGS=-O3 -Wall -I/usr/local/include -DHAVE_INLINE
LDFLAGS=-lgsl -lgslcblas -static
COMPILER=icc

evomod.exe: main.o evolution.o
	${COMPILER} -o evomod.exe main.o evolution.o $(LDFLAGS) $(CFLAGS)

burnin: main-burnin.o evolution.o
	${COMPILER} -o burnin.exe burnin.o evolution.o $(LDFLAGS) $(CFLAGS)

main.o: main.c evolution.o
	${COMPILER} -c main.c $(LDFLAGS) $(CFLAGS)

burnin.o: burnin.c evolution.o
	${COMPILER} -c burnin.c $(LDFLAGS) $(CFLAGS)

evolution.o : moments.o io.o
	${COMPILER} -c population/evolution.c $(LDFLAGS) $(CFLAGS)

io.o :
	${COMPILER} -c population/io.c $(LDFLAGS) $(CFLAGS)

moments.o :
	${COMPILER} -c population/moments.c $(LDFLAGS) $(CFLAGS)

clean:
	rm main.o io.o evolution.o moments.o pop.o

cleanall :
	rm main.o
