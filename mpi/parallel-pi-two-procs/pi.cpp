#include <cstdio>
#include <cmath>
#include <mpi.h>

constexpr int n = 840000;

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rank = 0;
  int nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  if(rank == 0)
    printf("Computing approximation to pi with N=%d\n", n);

  int step = n / nranks;
  int istart = 1 + rank * step;
  int istop = istart + step;

  double pi_local = 0.0;
  for (int i=istart; i <= istop; i++) {
    double x = (i - 0.5) / n;
    pi_local += 1.0 / (1.0 + x*x);
  }

  
  double totalpi = 0.0;
  MPI_Reduce(&pi_local, &totalpi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  totalpi *= 4.0 / n;

  if(rank == 0)
    printf("Approximate pi=%18.16f (exact pi=%10.8f)\n", totalpi, M_PI);
  MPI_Finalize();
}
