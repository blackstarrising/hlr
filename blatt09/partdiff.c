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

#include "partdiff.h"

struct calculation_arguments
{
  uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
  uint64_t  num_matrices;   /* number of matrices                             */
  double    h;              /* length of a space between two lines            */
  double    ***Matrix;      /* index matrix used for addressing M             */
  double    *M;             /* two matrices with real values                  */
  uint64_t  rank;
  uint64_t  world_size;
  uint64_t  Nh;             /* Heigth of a storage block                      */  
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
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
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;

	uint64_t N = arguments->N;
	uint64_t rank = arguments->rank;
	uint64_t world_size = arguments->world_size;

	uint64_t Nh_temp = (N-1)/world_size;
	if(rank < ((N-1) % world_size)){Nh_temp++;}
	arguments->Nh=Nh_temp;
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
	uint64_t const Nh = arguments->Nh;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (Nh + 2) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j < (Nh + 2); j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (Nh + 2)) + (j * (N + 1));
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
  uint64_t g, i, j, ih;                                /*  local variables for loops   */

	uint64_t const rank = arguments->rank;
	uint64_t const world_size = arguments->world_size;
	uint64_t const N = arguments->N;
	uint64_t const Nh = arguments->Nh;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
	  for (i = 0; i < (Nh + 2); i++)
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
	  int offset = ((N-1)/world_size*rank);
	  if((N-1)%world_size <= rank){offset += (N-1)%world_size;}
	  else if(rank < (N-1)%world_size){offset += rank;}
		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 0; i <= N; i++)
			{
			  if(rank == 0){Matrix[g][0][i] = 1.0 - (h * i);}
			  if(rank == world_size-1){Matrix[g][Nh+1][i] = h * i;}
			}
			for (ih = 0; ih <= Nh; ih++)
			{
			  Matrix[g][ih][0] = 1.0 - (h * (ih+offset));
			  Matrix[g][ih][N] = h * (ih+offset);
			}

			Matrix[g][Nh+1][0] = 0.0;
			Matrix[g][0][N] = 0.0;
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

/* **************************************************************************** */
/* calculateGaussSeidelMPI: solves the equation with Gauß-Seidel MPI-parallel   */
/* **************************************************************************** */
static
void
calculateGaussSeidelMPI (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
  uint64_t i, j;                                   /* local variables for loops */
  int m1 = 0;
  int m2 = 0;                                 /* used as indices for old and new matrices */
  double star;                                /* four times center value minus 4 neigh.b values */
  double residuum;                            /* residuum of current iteration */
  double maxresiduum;                         /* maximum residuum value of a slave in iteration */

  //Variables for parallelisation:
  short fertig = 0; //boolean to check finish status for precision

  uint64_t const rank = arguments->rank;
  uint64_t const world_size = arguments->world_size;
  uint64_t const Nh = arguments->Nh;

  uint64_t const N = arguments->N;
  double const h = arguments->h;

  double pih = 0.0;
  double fpisin = 0.0;

  int term_iteration = options->term_iteration;

  if (options->inf_func == FUNC_FPISIN)
    {
      pih = PI * h;
      fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

  if(rank == world_size-1){for(uint64_t i = 0; i < world_size; i++){MPI_Send(&fertig, 1, MPI_SHORT, 0, 0, MPI_COMM_WORLD);}}
  
  //TAGLIST
  //fertig: 0
  //Zeile von oben: 1
  //Maxresiduum: 2
  //Zeile von unten: 3
  while (term_iteration > 0)
    {
      if(rank == 0){MPI_Recv(&fertig, 1, MPI_SHORT, world_size-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);}
      else {MPI_Recv(&fertig, 1, MPI_SHORT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);}     
      
      /*if(fertig == 1){
	printf("Aborted die to reached precision limit");
	term_iteration = 0;
	break;
      }*/ // If the fertig is true, exit the while loop

      double** Matrix_u = arguments->Matrix[m1];
      
      maxresiduum = 0;

      if (rank == 0) //################################################################################
	{
	  /* over all rows */
	  for (i = 1; i <= Nh; i++)
	    {
	      double fpisin_i = 0.0;
	      
	      if (options->inf_func == FUNC_FPISIN)
		{
		  fpisin_i = fpisin * sin(pih * (double)i);
		}
	      
	      /* over all columns */
	      for (j = 1; j < N; j++)
		{
		  star = 0.25 * (Matrix_u[i-1][j] + Matrix_u[i][j-1] + Matrix_u[i][j+1] + Matrix_u[i+1][j]);
		  
		  if (options->inf_func == FUNC_FPISIN)
		    {
		      star += fpisin_i * sin(pih * (double)j);
		    }
		  
		  if (options->termination == TERM_PREC || term_iteration == 1)
		    {
		      residuum = Matrix_u[i][j] - star;
		      residuum = (residuum < 0) ? -residuum : residuum;
		      maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
		    }
		  
		  Matrix_u[i][j] = star;
		}
	    }

	  MPI_Send(&fertig, 1, MPI_SHORT, rank + 1, 0, MPI_COMM_WORLD);
	  MPI_Send(Matrix_u[Nh], N+1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
	  MPI_Send(&maxresiduum, 1, MPI_DOUBLE, rank + 1, 2, MPI_COMM_WORLD);
	  MPI_Recv(Matrix_u[Nh+1], N+1, MPI_DOUBLE, rank + 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
      else if (rank == world_size - 1) //###############################################################
	{
	  MPI_Recv(Matrix_u[0], N+1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&maxresiduum, 1, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  
	  //Only Calculate first row
	  i = 1;
	  double fpisin_i = 0.0;
	  
	  if (options->inf_func == FUNC_FPISIN)
	    {
	      fpisin_i = fpisin * sin(pih * (double)i);
	    }
	  
	  /* over all columns */
	  for (j = 1; j < N; j++)
	    {
	      star = 0.25 * (Matrix_u[i-1][j] + Matrix_u[i][j-1] + Matrix_u[i][j+1] + Matrix_u[i+1][j]);
	      
	      if (options->inf_func == FUNC_FPISIN)
		{
		  star += fpisin_i * sin(pih * (double)j);
		}
	      
	      if (options->termination == TERM_PREC || term_iteration == 1)
		{
		  residuum = Matrix_u[i][j] - star;
		  residuum = (residuum < 0) ? -residuum : residuum;
		  maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
		}
	      
	      Matrix_u[i][j] = star;
	    }
	  
	  MPI_Send(Matrix_u[1], N+1, MPI_DOUBLE, rank-1, 3, MPI_COMM_WORLD);
	  
	  //Iterate over rest of rows
	  /* over all rows */
	  for (i = 2; i <= Nh; i++)
	    {
	      if (options->inf_func == FUNC_FPISIN)
		{
		  fpisin_i = fpisin * sin(pih * (double)i);
		}
	      
	      /* over all columns */
	      for (j = 1; j < N; j++)
		{
		  star = 0.25 * (Matrix_u[i-1][j] + Matrix_u[i][j-1] + Matrix_u[i][j+1] + Matrix_u[i+1][j]);
		  
		  if (options->inf_func == FUNC_FPISIN)
		    {
		      star += fpisin_i * sin(pih * (double)j);
		    }
		  
		  if (options->termination == TERM_PREC || term_iteration == 1)
		    {
		      residuum = Matrix_u[i][j] - star;
		      residuum = (residuum < 0) ? -residuum : residuum;
		      maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
		    }
		  
		  Matrix_u[i][j] = star;
		}
	    }
	  
	  /* check for stopping calculation depending on termination method */
	  if (options->termination == TERM_PREC)
	    {
	      if (maxresiduum < options->term_precision)
		{
		  term_iteration = 0;
		  fertig = 1;
		}
	    }

	  MPI_Send(&fertig, 1, MPI_SHORT, 0, 0, MPI_COMM_WORLD);
	}
      
      else //###########################################################################################
	{
	  MPI_Recv(Matrix_u[0], N + 1, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(&maxresiduum, 1, MPI_DOUBLE, rank-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  //Only Calculate first row
	  i = 1;
	  double fpisin_i = 0.0;
	  
	  if (options->inf_func == FUNC_FPISIN)
	    {
	      fpisin_i = fpisin * sin(pih * (double)i);
	    }
	  
	  /* over all columns */
	  for (j = 1; j < N; j++)
	    {
	      star = 0.25 * (Matrix_u[i-1][j] + Matrix_u[i][j-1] + Matrix_u[i][j+1] + Matrix_u[i+1][j]);
	      
	      if (options->inf_func == FUNC_FPISIN)
		{
		  star += fpisin_i * sin(pih * (double)j);
		}
	      
	      if (options->termination == TERM_PREC || term_iteration == 1)
		{
		  residuum = Matrix_u[i][j] - star;
		  residuum = (residuum < 0) ? -residuum : residuum;
		  maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
		}
	      
	      Matrix_u[i][j] = star;
	    }

	  MPI_Send(Matrix_u[1], N+1, MPI_DOUBLE, rank-1, 3, MPI_COMM_WORLD);

		  
	  
	  //Iterate over rest of rows
	  /* over all rows */
	  for (i = 2; i <= Nh; i++)
	    {
	      if (options->inf_func == FUNC_FPISIN)
		{
		  fpisin_i = fpisin * sin(pih * (double)i);
		}
	      
	      /* over all columns */
	      for (j = 1; j < N; j++)
		{
		  star = 0.25 * (Matrix_u[i-1][j] + Matrix_u[i][j-1] + Matrix_u[i][j+1] + Matrix_u[i+1][j]);
		  
		  if (options->inf_func == FUNC_FPISIN)
		    {
		      star += fpisin_i * sin(pih * (double)j);
		    }
		  
		  if (options->termination == TERM_PREC || term_iteration == 1)
		    {
		      residuum = Matrix_u[i][j] - star;
		      residuum = (residuum < 0) ? -residuum : residuum;
		      maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
		    }
		  
		  Matrix_u[i][j] = star;
		}
	    }

	  MPI_Send(&fertig, 1, MPI_SHORT, rank + 1, 0, MPI_COMM_WORLD);
	  MPI_Send(Matrix_u[Nh], N+1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
	  MPI_Send(&maxresiduum, 1, MPI_DOUBLE, rank + 1, 2, MPI_COMM_WORLD);
	  MPI_Recv(Matrix_u[Nh+1], N+1, MPI_DOUBLE, rank + 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	} //####################################################################################################
      
      results->stat_iteration++;

      if (options->termination == TERM_ITER)
	{
	  term_iteration--;
	}

      //printf("Ich bin %li und laufe zum %i ten mal durch\n", rank, results->stat_iteration);
      
      
      if(fertig == 1){
	term_iteration = 0;
	if(rank == world_size - 1 && options->termination == TERM_PREC){
		for(uint64_t i = 0; i < world_size; i++){
			MPI_Send(Matrix_u[1], N+1, MPI_DOUBLE, rank-1, 3, MPI_COMM_WORLD); //Sende erneut die oberste Zeile an den Vorgänger, damit Prozesse beendet werden.
		}
	}
	break;
	} // If the fertig is true, exit the while loop

    }
  MPI_Bcast(&maxresiduum, 1, MPI_DOUBLE, world_size-1, MPI_COMM_WORLD);
  results->stat_precision = maxresiduum;
  
  results->m = m2;
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
DisplayMatrixLEGACY (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
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
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
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
  MPI_Init(&argc, &argv);
  int rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  struct options options;
  struct calculation_arguments arguments;
  struct calculation_results results;

  //Nur ein Prozess fragt optionen ab und teilt sie
  MPI_Datatype MPI_Struct_options;
  int structlen = 7;
  int blocklengths[structlen];
  MPI_Datatype types[structlen];
  MPI_Aint displacements[structlen];
  for (int i = 0; i < 6; i++)
    {
      blocklengths[i] = 1; types[i] = MPI_UNSIGNED_LONG;
    }
  blocklengths[6] = 1; types[6] = MPI_DOUBLE;

  displacements[0] = offsetof(options_struct, number);
  displacements[1] = offsetof(options_struct, method);
  displacements[2] = offsetof(options_struct, interlines);
  displacements[3] = offsetof(options_struct, inf_func);
  displacements[4] = offsetof(options_struct, termination);
  displacements[5] = offsetof(options_struct, term_iteration);
  displacements[6] = offsetof(options_struct, term_precision);

  MPI_Type_create_struct(structlen, blocklengths, displacements, types, &MPI_Struct_options);
  MPI_Type_commit(&MPI_Struct_options);

  if (rank == 0)
    {
      AskParams(&options, argc, argv);
    }
  MPI_Bcast(&options, 1, MPI_Struct_options, 0, MPI_COMM_WORLD);

  if(options.method == METH_JACOBI && world_size > 1)
    {
      if(rank == 0){printf("Aborted: Please use only one process for Jacobi Method!");}
      return 1;
    }
  
  arguments.rank = rank;
  arguments.world_size = world_size;
  initVariables(&arguments, &results, &options);  
  allocateMatrices(&arguments);
  
  MPI_Barrier(MPI_COMM_WORLD); //#####
  
  initMatrices(&arguments, &options);

  MPI_Barrier(MPI_COMM_WORLD); //#####

  gettimeofday(&start_time, NULL);
  if(world_size == 1){calculate(&arguments, &results, &options);}
  else {calculateGaussSeidelMPI(&arguments, &results, &options);}
  gettimeofday(&comp_time, NULL);

  MPI_Barrier(MPI_COMM_WORLD); //#####

  if(rank == 0){displayStatistics(&arguments, &results, &options);}

  MPI_Barrier(MPI_COMM_WORLD); //#####
 
  int N = arguments.N;
  int from = ((N-1)/world_size*rank) + 1;
  if((N-1)%world_size <= rank){from += (N-1)%world_size;}
  else if(rank < (N-1)%world_size){from += rank;}
  int to = from + arguments.Nh;
  
  if(world_size == 1){DisplayMatrixLEGACY(&arguments, &results, &options);}
  else{DisplayMatrix(&arguments, &results, &options, rank, world_size, from, to);}

  MPI_Barrier(MPI_COMM_WORLD); //#####
  
  freeMatrices(&arguments);

  MPI_Finalize();
  
  return 0;
}
