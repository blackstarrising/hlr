#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

#define MESSAGE_BUFFER_SIZE 50

int*
init (int N, int Nmod, int rank)
{
  //Every process has at least N numbers, some have N+1 (all in the beginning)
  //We allocate memory for n+1 for every process as every process can have n+1 numbers during runtime.
  int* buf = malloc(sizeof(int) * (N+1));

  srand(time(NULL)*(rank+1));//To make the different processes start at different point in the random number list, we multiply the "random" seed time by some factor individual to each process (rank)

  //All valid positions get filled random, all others filled with -1 to indicate that they are empty.
  for (int i = 0; i < N; i++)
    {
      // Do not modify "% 25"
      buf[i] = rand() % 25;
    }
  
  if((Nmod != 0) && (rank < Nmod))
    {
      buf[N] = rand() % 25;
    }
  else
    {
      buf[N] = -1;
    }
  return buf;
}

int*
circle (int* buf, int N, int rank, int world_size)
{
  //Send buf[0] of first process (rank 0) to last process (rank world_rank - 1)
  int endValue;
  if(rank == 0)
    {/*Send endvalue*/
      endValue=buf[0];
      MPI_Send(&endValue, 1, MPI_INT, world_size-1, 1, MPI_COMM_WORLD);
    }
  else if(rank == world_size-1)
    {/*recieve endvalue*/
      MPI_Recv(&endValue, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
  //Allocate memory to temporarily save two buffers (needed fpr exchanging)
  int* buf2 = malloc(sizeof(int) * (N+1));

  //Condition that specifies whether the terminating condition has been met
  int condition = 1;
  while(condition){
    for(int i = 0; i <= N; i++)
      {
	buf2[i] = buf[i]; //save "my" numbers
      }
    if(rank == 0) //Rank 0 sends first, then waits for recieve
      {
	MPI_Send(buf2, (N+1), MPI_INT, rank+1, 1, MPI_COMM_WORLD);
	MPI_Recv(buf, (N+1), MPI_INT, world_size-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    else if (rank == world_size - 1) //Finally, the last process recieves, then sends to 0, then checks the condition and sets it for himself if finished
      {
	MPI_Recv(buf, (N+1), MPI_INT, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Send(buf2, (N+1), MPI_INT, 0, 1, MPI_COMM_WORLD);
	if (buf[0] == endValue){condition = 0;}
      }
    else //All other ranks reviece, then send
      {
	MPI_Recv(buf, (N+1), MPI_INT, rank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Send(buf2, (N+1), MPI_INT, rank+1, 1, MPI_COMM_WORLD);
      }
    MPI_Bcast(&condition, 1, MPI_INT, world_size-1, MPI_COMM_WORLD); //Last Process shares his condition with all processes
  }
  free(buf2); //temp memory is freed
  return buf;
}

int
printConsole (int* buf, int N, int rank, int world_size)
{
  //Message is built from all entries in a buffer that are not empty (= -1)
  //To avoid appending to an empty buffer, this strange handling of the first entry (buf[0]) is needed. 
  char message[MESSAGE_BUFFER_SIZE*(N+1)];
  if(buf[0] != -1){sprintf(message, "rank %d: %d\n", rank, buf[0]);}
  else{sprintf(message, "\0");}
  for (int i = 1; i <= N; i++){
    if(buf[i] != -1){
      char buffer[MESSAGE_BUFFER_SIZE];
      sprintf(buffer, "rank %d: %d\n", rank, buf[i]);
      strcat(message, buffer);
    }
  }

  //All Processes send their output to process 0, process 0 prints all of them.
  //This whole problem could also be solved by just sending all numbers together with the corresponding ranks and just printing that. will probably change that up.
  if(rank == 0)
    {
      //Print own
      printf("%s", message);
      //Revieve and print others
      for(int i = 1; i < world_size; i++)
	{
	  MPI_Recv(message, MESSAGE_BUFFER_SIZE*(N+1), MPI_CHAR, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  printf("%s", message);
	}
    }
  else
    {
      MPI_Send(message, MESSAGE_BUFFER_SIZE*(N+1), MPI_CHAR, 0, 1, MPI_COMM_WORLD);
    }
  
  return EXIT_SUCCESS;
}


int
main (int argc, char** argv)
{
  if ((argc < 2))
    {
      printf("Arguments error!\nPlease specify a buffer size.\n");
      return EXIT_FAILURE;
    }
  //TODO all other strange cases have to be excluded (N = 0, N smaller than 0 etc...)

  //Init
  MPI_Init(&argc, &argv);
  
  int Ntot;
  int rank;
  int* buf;
  int world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  //total Array length
  Ntot = atoi(argv[1]);

  int Nmod = Ntot % world_size;
  int N = Ntot/world_size;
  
  buf = init(N, Nmod, rank);

  if(rank == 0){printf("\nBEFORE\n");}
  
  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################

  printConsole(buf, N, rank, world_size);

  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################
  
  circle(buf, N, rank, world_size);

  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################

  if(rank == 0){printf("\nAFTER\n");}
  
  printConsole(buf, N, rank, world_size);

  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################
  MPI_Finalize();

  free(buf);
  return EXIT_SUCCESS;
}
