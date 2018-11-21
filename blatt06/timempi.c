#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>


int main(int argc, char** argv)
{
  //Initialisiere MPI
  MPI_Init(NULL, NULL);

  //Lese Prozessrang aus
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  time_t rawtime;
  struct tm *time_info;
  char time_buffer[80];

  time(&rawtime);
  time_info = localtime(&rawtime);
  strftime(time_buffer, 80, "%x %X", time_info);

  struct timeval exacttime;
  gettimeofday(&exacttime, NULL);
  
  printf("%s: %s.%ld - rank: %d \n", processor_name, time_buffer, exacttime.tv_usec, world_rank);

  MPI_Finalize();
}
