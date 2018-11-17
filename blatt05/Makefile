# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff partdiff-posix
OBJS = partdiff.o partdiff-posix.o askparams.o

# Targets ...
all: partdiff partdiff-posix

partdiff: partdiff.o askparams.o Makefile
	$(CC) $(LFLAGS) -o $@ partdiff.o askparams.o $(LIBS)

partdiff-posix: partdiff-posix.o askparams.o Makefile
	$(CC) $(LFLAGS) -o $@ partdiff-posix.o askparams.o $(LIBS) -lpthread

partdiff.o: partdiff.c Makefile

partdiff-posix.o: partdiff-posix.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
