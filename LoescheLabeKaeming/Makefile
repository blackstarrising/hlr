# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff partdiff-openmp
OBJS = partdiff.o askparams.o

# Targets ...
all: partdiff partdiff-openmp

partdiff: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-openmp: partdiff-openmp.o askparams.o Makefile
	$(CC) $(LFLAGS) -fopenmp -o $@ partdiff-openmp.o askparams.o $(LIBS)
	
partdiff.o: partdiff.c
	$(CC) -c $(CFLAGS) $*.c

partdiff-openmp.o: partdiff-openmp.c 
	$(CC) -c $(CFLAGS) -fopenmp $*.c

askparams.o: askparams.c 
	$(CC) -c $(CFLAGS) $*.c
	
# Rule to create *.o from *.c
#%.o: %.c
#	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	rm -f partdiff-openmp.o
	$(RM) $(TGTS)
