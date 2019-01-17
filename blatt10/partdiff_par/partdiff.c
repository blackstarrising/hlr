/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and    **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <mpi.h>
#include<omp.h>

#include "partdiff.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

struct process_data
{
    int rank;
    int size;
    uint64_t chunk_size;
    uint64_t from;
    uint64_t to;
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options, struct process_data* process_data)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;

    //calculate data distibution
    uint64_t div = (arguments->N-1) / process_data->size;
    uint64_t rest = (arguments->N-1) % process_data->size;

    if(process_data->rank < rest)
    {
        process_data->chunk_size = div+3;
        process_data->from = process_data->rank * (div+1) + 1;
        process_data->to = process_data->from + div ;
    }
    else
    {
        process_data->chunk_size = div+2;
        process_data->from = process_data->rank * (div) + rest + 1;
        process_data->to = process_data->from + div - 1;
    }

    //printf("Rank %d, N-1: %d, from: %d, to: %d, chunk_size: %d\n",process_data->rank, arguments->N-1, process_data->from, process_data->to, process_data->chunk_size);

}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

static
void
allocateMatrices_mpi (struct calculation_arguments* arguments, struct process_data* process_data)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;
    uint64_t const chunk_size = process_data->chunk_size;

	arguments->M = allocateMemory(arguments->num_matrices * chunk_size * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory(chunk_size * sizeof(double*));

		for (j = 0; j < chunk_size; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * chunk_size * (N + 1)) + (j * (N + 1));
		}
	}
}


/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t g, i, j;                                /*  local variables for loops   */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 0; i <= N; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				Matrix[g][0][i] = 1.0 - (h * i);
				Matrix[g][N][i] = h * i;
			}

			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}

static
void
initMatrices_mpi (struct calculation_arguments* arguments, struct options const* options, struct process_data* process_data)
{
	uint64_t g, i, j;                                /*  local variables for loops   */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

    uint64_t chunk_size = process_data->chunk_size;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i < chunk_size; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{

            for(i = 0; i < chunk_size; i++)
            {
				Matrix[g][i][0] = 1.0 - (h * (process_data->from + i - 1));
				Matrix[g][i][N] = h * (process_data->from + i - 1);
            }
            
			for (i = 0; i <= N; i++)
			{
                if(process_data->rank == 0)
                {
                    Matrix[g][0][i] = 1.0 - (h * i);
                }
                else if(process_data->rank == process_data->size-1)
                {
                    Matrix[g][chunk_size-1][i] = h * i;
                }
			}

            if(process_data->rank == 0)
            {
                Matrix[g][0][N] = 0.0;
            }
            else if(process_data->rank == process_data->size-1)
            {
                Matrix[g][chunk_size-1][0] = 0.0;
            }
		}
	}
}


/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;                                   /* local variables for loops */
	int m1, m2;                                 /* used as indices for old and new matrices */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxresiduum = 0;

		/* over all rows */
		for (i = 1; i < N; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxresiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;
}

static
void
calculate_jacobi_mpi (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct process_data* process_data)
{
	int i, j;                                   /* local variables for loops */
	int m1, m2;                                 /* used as indices for old and new matrices */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

    m1 = 0;
    m2 = 1;

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}


    int rank = process_data->rank;
    int size = process_data->size;

    int proc_next = rank + 1;
    int proc_before = rank - 1;

    int from = process_data->from;
    int chunk_size = process_data->chunk_size;

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxresiduum = 0;

        MPI_Request req;

        if(proc_next <= size-1)
        {
            MPI_Isend(&Matrix_In[chunk_size-2][0], N+1, MPI_DOUBLE, proc_next, 0, MPI_COMM_WORLD, &req);
            MPI_Recv(&Matrix_In[chunk_size-1][0], N+1, MPI_DOUBLE, proc_next, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    
        if(proc_before >= 0 )
        {
            MPI_Isend(&Matrix_In[1][0], N+1, MPI_DOUBLE, proc_before, 0, MPI_COMM_WORLD, &req);
            MPI_Recv(&Matrix_In[0][0], N+1, MPI_DOUBLE, proc_before, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);

		/* over all rows */
		for (i = 1; i < chunk_size-1 ; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)(i + from - 1));
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		results->stat_iteration++;

        double tmp;
        MPI_Allreduce(&maxresiduum, &tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        maxresiduum = tmp;

		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxresiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;
}

static
void
calculate_gauss_seidel_mpi (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct process_data* process_data)
{
	int i, j;                                   /* local variables for loops */
	int m1, m2;                                 /* used as indices for old and new matrices */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

    m1 = 0;
    m2 = 0;

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

    int rank = process_data->rank;
    int size = process_data->size;

    int proc_next = rank + 1;
    int proc_before = rank - 1;

    int from = process_data->from;
    int chunk_size = process_data->chunk_size;

    MPI_Request req;

    int stop = 0;
    int detect_stop = 0;

    double** Matrix_Out = arguments->Matrix[m1];
    double** Matrix_In  = arguments->Matrix[m2];
    if(proc_before >= 0 )
    {
        MPI_Isend(&Matrix_In[1][0], N+1, MPI_DOUBLE, proc_before, 0, MPI_COMM_WORLD, &req);
    }

	while (term_iteration > 0 && stop != 1)
	{
		maxresiduum = 0;

		/* over all rows */
		for (i = 1; i < chunk_size-1 ; i++)
		{
            if (i==1)
            {
                if(proc_before >= 0 )
                {
                    MPI_Recv(&Matrix_In[0][0], N+1, MPI_DOUBLE, proc_before, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&maxresiduum, 1, MPI_DOUBLE, proc_before, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if (options->termination == TERM_PREC)
                    {
                        MPI_Recv(&stop, 1, MPI_INT, proc_before, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
                else
                {
                    if (options->termination == TERM_PREC && results->stat_iteration >= size)
                    {
                        MPI_Recv(&stop, 1, MPI_INT, size-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
            else if (i == chunk_size-2 && proc_next <= size-1)
            {
                MPI_Recv(&Matrix_In[chunk_size-1][0], N+1, MPI_DOUBLE, proc_next, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)(i + from - 1));
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}

            if (i==1)
            {
                if(proc_before >= 0 )
                {
                    MPI_Isend(&Matrix_In[1][0], N+1, MPI_DOUBLE, proc_before, 0, MPI_COMM_WORLD, &req);
                }
            }
            else if (i == chunk_size-2)
            {
                if(proc_next <= size-1)
                {
                    MPI_Isend(&Matrix_In[chunk_size-2][0], N+1, MPI_DOUBLE, proc_next, 0, MPI_COMM_WORLD, &req);
                    MPI_Isend(&maxresiduum, 1, MPI_DOUBLE, proc_next, 0, MPI_COMM_WORLD, &req);
                    if (options->termination == TERM_PREC)
                    {
                        MPI_Isend(&stop, 1, MPI_INT, proc_next, 0, MPI_COMM_WORLD, &req);
                    }
                }
            }
		}

		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
            if(detect_stop == 0 && rank == size-1)
            {
                if (maxresiduum < options->term_precision && rank == size-1)
                {
                    detect_stop = 1;
                }
                MPI_Isend(&detect_stop, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &req);
            }
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;
    MPI_Bcast(&results->stat_precision, 1, MPI_DOUBLE, size-1, MPI_COMM_WORLD);
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion DisplayMatrix:                               **/
/**                                                                        **/
/** Die Funktion DisplayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf ("\n");
	}

	fflush (stdout);
}

/**
 * rank and size are the MPI rank and size, respectively.
 * from and to denote the global(!) range of lines that this process is responsible for.
 *
 * Example with 9 matrix lines and 4 processes:
 * - rank 0 is responsible for 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7.
 *   Lines 0 and 8 are not included because they are not calculated.
 * - Each process stores two halo lines in its matrix (except for ranks 0 and 3 that only store one).
 * - For instance: Rank 2 has four lines 0-3 but only calculates 1-2 because 0 and 3 are halo lines for other processes. It is responsible for (global) lines 5-6.
 */
static
void
DisplayMatrix_mpi (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
  int const elements = 8 * options->interlines + 9;

  int x, y;
  double** Matrix = arguments->Matrix[results->m];
  MPI_Status status;

  /* first line belongs to rank 0 */
  if (rank == 0)
    from--;

  /* last line belongs to rank size - 1 */
  if (rank + 1 == size)
    to++;

  if (rank == 0)
    printf("Matrix:\n");

  for (y = 0; y < 9; y++)
  {
    int line = y * (options->interlines + 1);

    if (rank == 0)
    {
      /* check whether this line belongs to rank 0 */
      if (line < from || line > to)
      {
        /* use the tag to receive the lines in the correct order
         * the line is stored in Matrix[0], because we do not need it anymore */
        MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
      }
    }
    else
    {
      if (line >= from && line <= to)
      {
        /* if the line belongs to this process, send it to rank 0
         * (line - from + 1) is used to calculate the correct local address */
        MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
      }
    }

    if (rank == 0)
    {
      for (x = 0; x < 9; x++)
      {
        int col = x * (options->interlines + 1);

        if (line >= from && line <= to)
        {
          /* this line belongs to rank 0 */
          printf("%7.4f", Matrix[line][col]);
        }
        else
        {
          /* this line belongs to another rank and was received above */
          printf("%7.4f", Matrix[0][col]);
        }
      }

      printf("\n");
    }
  }

  fflush(stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{

	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;
    struct process_data process_data;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_data.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &process_data.size);

    if(process_data.rank == 0)
    {
        AskParams(&options, argc, argv);
    }

    MPI_Bcast(&options, sizeof(struct options), MPI_BYTE, 0, MPI_COMM_WORLD);

	initVariables(&arguments, &results, &options, &process_data);

    if(process_data.size == 1 || (process_data.size > arguments.N-1 ))
    {
        if(process_data.rank == 0)
        {
            allocateMatrices(&arguments);
            initMatrices(&arguments, &options);

            gettimeofday(&start_time, NULL);
            calculate(&arguments, &results, &options);
            gettimeofday(&comp_time, NULL);

            displayStatistics(&arguments, &results, &options);
            DisplayMatrix(&arguments, &results, &options);

            freeMatrices(&arguments);
        }
    }
    else if (options.method == METH_JACOBI)
    {
        allocateMatrices_mpi(&arguments, &process_data);
        initMatrices_mpi(&arguments, &options, &process_data);

        gettimeofday(&start_time, NULL);
        calculate_jacobi_mpi(&arguments, &results, &options, &process_data);
        gettimeofday(&comp_time, NULL);

        if(process_data.rank == 0)
        {
            displayStatistics(&arguments, &results, &options);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        DisplayMatrix_mpi(&arguments, &results, &options, process_data.rank, process_data.size, process_data.from, process_data.to);

        freeMatrices(&arguments);
    }
    else if (options.method == METH_GAUSS_SEIDEL)
    {
        allocateMatrices_mpi(&arguments, &process_data);
        initMatrices_mpi(&arguments, &options, &process_data);

        gettimeofday(&start_time, NULL);
        calculate_gauss_seidel_mpi(&arguments, &results, &options, &process_data);
        gettimeofday(&comp_time, NULL);

        if(process_data.rank == 0)
        {
            displayStatistics(&arguments, &results, &options);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        DisplayMatrix_mpi(&arguments, &results, &options, process_data.rank, process_data.size, process_data.from, process_data.to);

        freeMatrices(&arguments);
    }


    MPI_Finalize();

	return 0;
}
