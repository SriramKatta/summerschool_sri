#include <cstdio>
#include <mpi.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    int rank, thread_id;
    int provided, required=MPI_THREAD_FUNNELED;

    MPI_Init_thread(&argc, &argv, required, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#pragma omp parallel private(thread_id)
    {
        thread_id = omp_get_thread_num();
        printf("I'm thread %d in process %d\n", thread_id, rank);
    }

    if (rank == 0) {
        printf("\nProvided thread support level: %d\n", provided);
        printf("  %d - MPI_THREAD_SINGLE\n", MPI_THREAD_SINGLE);
        printf("  %d - MPI_THREAD_FUNNELED\n", MPI_THREAD_FUNNELED);
        printf("  %d - MPI_THREAD_SERIALIZED\n", MPI_THREAD_SERIALIZED);
        printf("  %d - MPI_THREAD_MULTIPLE\n", MPI_THREAD_MULTIPLE);
    }
    MPI_Finalize();
    return 0;
}
