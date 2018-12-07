# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff-par partdiff-par-hybrid
OBJS1 = partdiff.o askparams.o
OBJS2 = partdiff-hybrid.o askparams.o

# Targets ...
all: partdiff-par partdiff-par-hybrid

partdiff-par: $(OBJS1) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS1) $(LIBS)

partdiff-par-hybrid: $(OBJS2) Makefile
	$(CC) $(LFLAGS) -fopenmp -o $@ $(OBJS2) $(LIBS)

partdiff.o: partdiff.c
	$(CC) -c $(CFLAGS) partdiff.c

partdiff-hybrid.o: partdiff.c
	$(CC) -c $(CFLAGS) -fopenmp -o partdiff-hybrid.o partdiff.c

askparams.o: askparams.c
	$(CC) -c $(CFLAGS) askparams.c

clean:
	$(RM) $(OBJS1)
	$(RM) $(OBJS2)
	$(RM) $(TGTS)
