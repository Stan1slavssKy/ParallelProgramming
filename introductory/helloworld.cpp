#include <mpi.h>
#include <stdio.h>
#include <sched.h>

int main( int argc , char **argv )
{
    MPI_Init(NULL, NULL);

    int commsize = 0;
    int task_id = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);

    printf("Hello world! [task_id = %d out of %d processors, run on cpu = %d]\n", task_id, commsize, sched_getcpu());

    MPI_Finalize();
}