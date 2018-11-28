#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

#define MESSAGE_BUFFER_SIZE 50

int*
init (int N, int Nlast, int last)
{
	// TODO
	int* buf = malloc(sizeof(int) * N);

	srand(time(NULL)*last);

	if((Nlast != 0) && (last == 1))
	  {
	    for (int i = 0; i < Nlast; i++)
	      {
		// Do not modify "% 25"
		buf[i] = rand() % 25;
	      }
	    for (int i = Nlast; i < N; i++)
	      {
		buf[i] = -1;
	      }
	  }
	else
	  {
	    for (int i = 0; i < N; i++)
	      {
		// Do not modify "% 25"
		buf[i] = rand() % 25;
	      }
	  }
	return buf;
}

int*
circle (int* buf)
{
	// TODO
	return buf;
}

int
main (int argc, char** argv)
{
  if ((argc < 2))
    {
      printf("Arguments error!\nPlease specify a buffer size.\n");
      return EXIT_FAILURE;
    }

  //Init
  MPI_Init(&argc, &argv);
  
  int Ntot;
  int rank;
  int* buf;
  int world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if ((argc < 2))
    {
      if(rank == 0){printf("Arguments error!\nPlease specify a buffer size.\n");}
      return EXIT_FAILURE;
    }

  
  //total Array length
  Ntot = atoi(argv[1]);

  int Nlast = Ntot % world_size;
  int N = Ntot/world_size;
  
  buf = init(N, Nlast, world_size - rank);

  if(rank == 0){printf("\nBEFORE\n");}
  
  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################

  char message[MESSAGE_BUFFER_SIZE*N];
  sprintf(message, "rank %d: %d\n", rank, buf[0]);
  for (int i = 1; i < N; i++){
    char buffer[MESSAGE_BUFFER_SIZE];
    sprintf(buffer, "rank %d: %d\n", rank, buf[i]);
    strcat(message, buffer);
  }
  
  if(rank == 0)
    {
      //Print own
      printf("%s", message);
      //Revieve and print others
      for(int i = 1; i < world_size; i++)
	{
	  MPI_Recv(message, MESSAGE_BUFFER_SIZE*N, MPI_CHAR, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  printf("%s", message);
	}
    }
  else
    {
      MPI_Send(message, MESSAGE_BUFFER_SIZE*N, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
    }

  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################
  
  circle(buf);

  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################
  
  if(rank == 0){printf("\nAFTER\n");}
  sprintf(message, "rank %d: %d\n", rank, buf[0]);
  for (int i = 1; i < N; i++){
    char buffer[MESSAGE_BUFFER_SIZE];
    sprintf(buffer, "rank %d: %d\n", rank, buf[i]);
    strcat(message, buffer);
  }
  
  if(rank == 0)
    {
      //Print own
      printf("%s", message);
      //Revieve and print others
      for(int i = 1; i < world_size; i++)
	{
	  MPI_Recv(message, MESSAGE_BUFFER_SIZE*N, MPI_CHAR, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  printf("%s", message);
	}
    }
  else
    {
      MPI_Send(message, MESSAGE_BUFFER_SIZE*N, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
    }

  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################
  MPI_Finalize();

  free(buf);
  return EXIT_SUCCESS;
}
