EXECS=timempi
MPICC?=mpicc

all: ${EXECS}

timempi: timempi.c
	${MPICC} -o timempi timempi.c

clean:
	rm ${EXECS}
