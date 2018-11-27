#include <stdio.h>
#include <mpi.h>
#include <limits.h>
//Headerfile für gettimeofday
#include <sys/time.h>
#include <time.h>
//Festlegen einiger Konstanten
#define TIME_BUFFER_SIZE 50
#define NODE_MESSAGE_BUFFER_SIZE 100
#define MASTER_RANK 0
#define NODE_MESSAGE_TAG 0

int main(int argc, char **argv)
{
	int ierr, rank, world_size;
	ierr = MPI_Init(&argc, &argv);

	//Rang des Prozess auslesen.
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Vorberiten der Variable für die Nachricht von den Nodes
	char node_message[NODE_MESSAGE_BUFFER_SIZE];
  // usecs wird hier von jedem Prozess gebraucht
  long usecs, lowest_usecs;
	if (rank == MASTER_RANK) {
		for(int i = 1; i < world_size; i++) {
			//MPI Nachrichten empfangen
			MPI_Recv(node_message, NODE_MESSAGE_BUFFER_SIZE, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("%s", node_message);
		}
    usecs = LONG_MAX;
	} else {
		//Inititalisieren der Variable für den Hostname
		char hostname[MPI_MAX_PROCESSOR_NAME];
		int processor_name_length;
		MPI_Get_processor_name(hostname, &processor_name_length);

		// Variablen für die Zeit
		struct timeval current_time;
		struct tm *local_time;
		char time_chars[TIME_BUFFER_SIZE];

		gettimeofday(&current_time,(struct timezone *)0);
		local_time = localtime(&current_time.tv_sec);
		usecs = current_time.tv_usec;
		strftime(time_chars, sizeof(time_chars), "%F %T", local_time);

		//Nachricht zusammensetzen
		sprintf(node_message, "%s:\t %s.%i\n", hostname, time_chars, usecs);
		//Senden an den Master
		MPI_Send(node_message, NODE_MESSAGE_BUFFER_SIZE, MPI_CHAR, MASTER_RANK, NODE_MESSAGE_TAG, MPI_COMM_WORLD);
	}

  //Minimum von usec
  MPI_Reduce(&usecs, &lowest_usecs, 1, MPI_LONG, MPI_MIN, MASTER_RANK, MPI_COMM_WORLD);

  if(rank == MASTER_RANK) {
    printf("Kleinste Millisekunden: %i\n", lowest_usecs);
  }
	//Auf alle Prozesse warten
	MPI_Barrier(MPI_COMM_WORLD);

	printf("Rang %i beendet jetzt!\n", rank);

	ierr = MPI_Finalize();

	return ierr;
}
