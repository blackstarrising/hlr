MPICC?=mpicc

all: timempi timempi2

timempi: timempi.c
	${MPICC} -g -o timempi timempi.c

timempi2: timempi2.c
	${MPICC} -g -o timempi2 timempi2.c

clean:
	rm timempi
	rm timempi2
