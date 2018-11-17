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
#define NTHREADS 12
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <pthread.h>
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

// @EDIT Dieses Struct enthält alle relevanten Daten für den Thread
struct step_params
{
	uint64_t 	threadid;
	double		**input;
	double		**output;
	uint64_t 	pars_term_iter;
	uint64_t 	pars_inf_func;
	uint64_t 	pars_termination;

};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */
// @EDIT Für unsere Thread-Funktion müssen wir einige zusätzliche globale
// variablen definieren
// @EDIT Anzahl Iterationen, die EIN thread machen muss:
int num_iterations;
// @EDIT: Da nicht alle Threads immer perfekt aufgeteilt werden können,
// müssen einige evtl eine Iteration mehr machen:
int num_rest;
// @EDIT Anzahl Elemente pro Zeile:
int num_elements;
// @EDIT Das ist ein Array, wo jeder Thread das maximale residuum der von ihm
// bearbeiteten Elemente hereinschreibt (damit sparen wir uns ein mutex):
static double maxres[NTHREADS];
double pih;
double fpisin;

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


void *
calculate_step(void* pars)
{
	// @EDIT Zunächst legen wir einige Startparameter fest
	struct step_params* mypars;
	mypars = (struct step_params*) pars;
	int my_id = mypars->threadid;
	double** my_input = mypars->input;
	double** my_output = mypars->output;
	double my_res = 0.0;
	double my_maxres = 0.0;
	int start_i = my_id;
	int start_j = 1;
	int my_iter = num_iterations;
	double star;
	//@EDIT Definiere ob dieser Thread eine Iteration mehr machen muss
	if(my_id <= num_rest)
	{
		my_iter += 1;
	}
	// @EDIT Hier wird der Fall abgefangen, dass es mehr threads geben kann,
	// als Elemente in der Zeile sind. Dann muss das i angepasst werden und
	// j muss entsprechend erhöht werden.
	while (num_elements < start_i)
	{
		start_i = start_i - num_elements;
		start_j += 1;
	}
	// @EDIT Das ist nun die kopierte for-loop aus dem Sequenziellen Programm;
	// natürlich sind einige Parameternamen angepasst worden.
	for (int count = 0; count < my_iter; count++)
	{
		double fpisin_i = 0.0;

		star = 0.25 * (my_input[start_i-1][start_j] + my_input[start_i][start_j-1] + my_input[start_i][start_j+1] + my_input[start_i+1][start_j]);

		if (mypars->pars_inf_func == FUNC_FPISIN)
		{
			fpisin_i = fpisin * sin(pih * (double)start_i);
			star += fpisin_i * sin(pih * (double)start_j);
		}

		if (mypars->pars_termination == TERM_PREC || mypars->pars_term_iter == 1)
		{
			my_res = my_input[start_i][start_j] - star;
			my_res = (my_res < 0) ? -my_res : my_res;
			my_maxres = (my_res < my_maxres) ? my_maxres : my_res;
		}

		my_output[start_i][start_j] = star;
		// @EDIT Wenn dieses Element berechnet ist, geht der thread mit der Schrittgröße
		// der Anzahl an threads weiter durch die Matrix. Sollte der Fall eintreten, dass
		// der Thread die Zeile überschreitet, werden die i und j-Werte wieder angepasst (wie oben)
		start_i += NTHREADS;
		while(num_elements < start_i)
		{
			start_i = start_i - num_elements;
			start_j += 1;
		}
	}
	// @EDIT Das maximale Residuum dieses threads wird in den entspechenden Teil
	// des globalen maxres arrays geschrieben
	maxres[my_id-1] = my_maxres;
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	// @EDIT p brauchen wir später nur um die Matrizen zu tauschen.
	int p;
	int m1, m2;                                 /* used as indices for old and new matrices */
	static double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	pih = 0.0;
	fpisin = 0.0;

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

		maxresiduum = 0.0;
		// @EDIT Hier definieren wir unsere pthreads und erzeugen außerdem ein
		// joinable Attribut.
		pthread_t threads[NTHREADS];
		pthread_attr_t attr;
		pthread_attr_init(&attr);
   	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		void* status;
		// @EDIT In diesem Array sind die structs mit den Daten für die Threads
		struct step_params data_array[NTHREADS];

		// @EDIT Sollten die Interlines nicht durch die Thread-Zahl teilbar sein,
		// müssen später einige threads ein Element mehr als andere Threads
		// ausrechnen
		num_rest = ((N-1)*(N-1))%NTHREADS;
		// @EDIT So viele Iterationen macht ein Thread (plus evtl eine mehr)
		num_iterations = (((N-1)*(N-1))-num_rest)/NTHREADS;
		num_elements = N-1;

		for (int k = 0; k < NTHREADS; k++)
		{
			// @EDIT Hier legen wir für jeden Struct die Daten für den Thread fest
			(&data_array[k])->threadid = k+1;
			(&data_array[k])->input = Matrix_In;
			(&data_array[k])->output = Matrix_Out;
			(&data_array[k])->pars_term_iter = term_iteration;
			(&data_array[k])->pars_termination = options->termination;
			(&data_array[k])->pars_inf_func = options->inf_func;
			// @EDIT Nun erzeugen wir die Threads mit den Daten und rufen damit die
			// Funktion "calculate_step" auf
			pthread_create(&threads[k], &attr, calculate_step,(void *) &data_array[k]);
		}

		// @EDIT Jetzt wird das Attribut "zerstört" und die Threads werden wieder
		// gejoined.
		pthread_attr_destroy(&attr);
 	  for(int k = 0; k < NTHREADS; k++)
		{
			pthread_join(threads[k], &status);
		}
		// @EDIT Aus dem "maxres" Array wird nun noch das tatsächliche maximale
		// Residuum ermittelt.
		for (int l = 0; l < NTHREADS; l++)
		{
			maxresiduum = (maxres[l] > maxresiduum) ? maxres[l] : maxresiduum;
		}

		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		p = m1;
		m1 = m2;
		m2 = p;

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

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	AskParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
	calculate(&arguments, &results, &options);
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	DisplayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);
	// @EDIT Das soll man laut wiki am Ende der main machen ;)
	pthread_exit(NULL);

	return 0;
}
