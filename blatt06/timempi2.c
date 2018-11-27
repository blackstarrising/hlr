#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
  //Initialisiere MPI
  MPI_Init(NULL, NULL);

  //Lese Anzahl Prozesse aus
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  //Lese Prozessrang aus
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  char outputstring[50];
  int usecs = 1000000;//Da Rang 0 diesen WErt nicht ändert muss er dort (durch diese initialisierung) größer sein als alle möglichen werte der anderen Prozesse.

  MPI_Barrier(MPI_COMM_WORLD);

  if(world_rank == 0)
    {
      //Prozess 0 empfängt die Strings aller anderen Prozesse.
      for (int j = 1; j < world_size; j++)
	{
	  MPI_Recv(outputstring, 50, MPI_CHAR, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  printf("%s\n", outputstring);
	}
    }
  else
    {
      //Generiere Prozessorname
      char processor_name[MPI_MAX_PROCESSOR_NAME];
      int namelen;
      MPI_Get_processor_name(processor_name, &namelen);
      
      //Generiere Timestamp
      time_t rawtime;
      struct tm *time_info;
      char time_buffer[80];
      time(&rawtime);
      time_info = localtime(&rawtime);
      strftime(time_buffer, 30, "%x %X", time_info);
      struct timeval exacttime;
      gettimeofday(&exacttime, NULL);
      usecs = (int) exacttime.tv_usec;
      
      //Generiere Outputstring
      sprintf(outputstring, "%s: %s.%i", processor_name, time_buffer, usecs);
      //sprintf(outputstring,"test");
      
      //Sende Outputstring an Prozess 0
      MPI_Send(outputstring, 50, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
  
  MPI_Barrier(MPI_COMM_WORLD); //#######################################################
  
  //Druck der kleinsten usecs mit Reduce:
  int min_usecs;
  MPI_Reduce(&usecs, &min_usecs, world_size-1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
  if(world_rank == 0){printf("%i\n", min_usecs);}
  
  MPI_Barrier(MPI_COMM_WORLD); //#######################################################
  
  printf("Rang %i beendet jetzt!\n", world_rank);
  
  MPI_Finalize();
}
