// Main solver routines for heat equation solver

#include <mpi.h>

#include "heat.hpp"

// Exchange the boundary values
void exchange(Field &field, const ParallelData parallel)
{
  MPI_Request req[4];
  MPI_Status status[4];
  int top = 0;
  int bot = field.nx + 1;

  int dombot = bot - 1;
  int domtop = top + 1;

  // TODO start: implement halo exchange

  // You can utilize the data() method of the Matrix class to obtain pointer
  // to element, e.g. field.temperature.data(i, j)
  // Send to up, receive from down
  double *sbuf = field.temperature.data(domtop, 0);
  double *rbuf = field.temperature.data(bot, 0);
    MPI_Isend(sbuf, field.ny + 2, MPI_DOUBLE, parallel.nup, 0, MPI_COMM_WORLD, &req[0]);
    MPI_Irecv(rbuf, field.ny + 2, MPI_DOUBLE, parallel.ndown, 0, MPI_COMM_WORLD, &req[1]);
  
  // Send to down, receive from up
  double *sbuf2 = field.temperature.data(dombot, 0);
  double *rbuf2 = field.temperature.data(top, 0);
  MPI_Isend(sbuf2, field.ny + 2, MPI_DOUBLE, parallel.ndown, 1, MPI_COMM_WORLD, &req[2]);
  MPI_Irecv(rbuf2, field.ny + 2, MPI_DOUBLE, parallel.nup, 1, MPI_COMM_WORLD, &req[3]);
  
  MPI_Waitall(4,req, status);

  // TODO end
}

// Update the temperature values using five-point stencil */
void evolve(Field &curr, const Field &prev, const double a, const double dt)
{

  // Compilers do not necessarily optimize division to multiplication, so make it explicit
  auto inv_dx2 = 1.0 / (prev.dx * prev.dx);
  auto inv_dy2 = 1.0 / (prev.dy * prev.dy);

  // Determine the temperature field at next time step
  // As we have fixed boundary conditions, the outermost gridpoints
  // are not updated.
  for (int i = 1; i < curr.nx + 1; i++)
  {
    for (int j = 1; j < curr.ny + 1; j++)
    {
      curr(i, j) = prev(i, j) + a * dt * ((prev(i + 1, j) - 2.0 * prev(i, j) + prev(i - 1, j)) * inv_dx2 + (prev(i, j + 1) - 2.0 * prev(i, j) + prev(i, j - 1)) * inv_dy2);
    }
  }
}
