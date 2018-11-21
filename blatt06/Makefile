EXECS=timempi
MPICC?=mpicc

all: ${EXECS}

timempi: timempi.c
	${MPICC} -g -o timempi timempi.c

clean:
	rm ${EXECS}
