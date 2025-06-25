#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    int nranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
    char hostname[500];
    int hostnamelen = 0;
    MPI_Get_processor_name(hostname, &hostnamelen);
    std::cout << "Hello from rank " << rank << " with name " << hostname << std::endl;

    if(0 ==  rank){
        std::cout << "total number of mpi processes " << nranks << std::endl;
    }
    
    if(rank == (nranks-1)){
        std::cout << "I am last but not least" << std::endl;
    }

    if(rank == 42){
        std::cout << "I'm the Answer to the Ultimate Question of Life, the Universe, and Everything!" << std::endl;
    }
    MPI_Finalize();
}
