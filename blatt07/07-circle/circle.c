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

//Two different print versions!
//Version 1: Each Process generates a String, sends it to 0 who prints it.
//Pro: 0 never has all array numbers
//Con: Char array send is not as great as just sending an int array.
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

//Version 2: Each process just sends their buffer to 0, 0 prints them
//Pro: shorter, no strange string appending required, easier send
//Con: Process 0 has all numbers from the array (although not at the same time!)
//Wr prefer version 2 but left version 1 in here if V2 is violating the rule to not have all numbers known to one process
int
printConsole2 (int* buf, int N, int rank, int world_size)
{
  int* bufp = malloc(sizeof(int) * (N+1));//Need another temp array to not overwrite buf (as it is used afterwards)
  if(rank == 0)
    {
      //Print own
      for (int i = 0; i <= N; i++)
	{
	  if(buf[i] != -1){printf("rank 0: %d\n",buf[i]);}
	}
      //Recieve and print others
      for(int i = 1; i < world_size; i++)
	{
	  MPI_Recv(bufp, (N+1), MPI_INT, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  for (int j = 0; j <= N; j++)
	    {
	      if(bufp[j] != -1){printf("rank %d: %d\n", i, bufp[j]);}
	    }
	}
    }
  else
    {
      MPI_Send(buf, (N+1), MPI_INT, 0, 2, MPI_COMM_WORLD);
    }
  free(bufp);
  return EXIT_SUCCESS;
}

int
main (int argc, char** argv)
{
  //Init
  MPI_Init(&argc, &argv);
  
  int Ntot;
  int rank;
  int* buf;
  int world_size;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  //TODO catch nun-number inputs
  if ((argc < 2))
    {
      if(rank == 0){printf("Arguments error!\nPlease specify a valid buffer size.\n");}
      return EXIT_FAILURE;
    }

  //total Array length
  Ntot = atoi(argv[1]);
  
  if (Ntot < 1)
    {
      if(rank == 0){printf("Arguments error!\nPlease specify a buffer size > 0.\n");}
      return EXIT_FAILURE;
    }
  
  int Nmod = Ntot % world_size;
  int N = Ntot/world_size;
  
  buf = init(N, Nmod, rank);

  if(rank == 0){printf("\nBEFORE\n");}
  
  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################
  
  printConsole2(buf, N, rank, world_size);

  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################
  
  circle(buf, N, rank, world_size);

  if(rank == 0){printf("\nAFTER\n");}

  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################
  
  printConsole2(buf, N, rank, world_size);
  
  MPI_Barrier(MPI_COMM_WORLD); //#####################################################################
  MPI_Finalize();
  
  free(buf);
  return EXIT_SUCCESS;
}
