/* I/O related functions for heat equation solver */

#include <string>
#include <iomanip>
#include <fstream>
#include <string>
#include <mpi.h>

#include "matrix.hpp"
#include "heat.hpp"
#include "pngwriter.h"

// Write a picture of the temperature field
void write_field(const Field& field, const int iter, const ParallelData parallel)
{

    auto height = field.nx * parallel.size;
    auto width = field.ny;
    auto full_data = Matrix<double>(height, width);

    // array for MPI sends and receives
    auto tmp_mat = Matrix<double> (field.nx, field.ny);
    // collect data
    for (int i = 0; i < field.nx; i++)
        for (int j = 0; j < field.ny; j++)
                tmp_mat(i, j) = field(i + 1, j + 1);

    int numelems = field.nx * field.ny;
    
    MPI_Gather(tmp_mat.data(), numelems, MPI_DOUBLE, full_data.data(), numelems, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (0 == parallel.rank){
        // Write out the data to a png file
        std::ostringstream filename_stream;
        filename_stream << "heat_" << std::setw(4) << std::setfill('0') << iter << ".png";
        std::string filename = filename_stream.str();
        save_png(full_data.data(), height, width, filename.c_str(), 'c');
    }


}

// Read the initial temperature distribution from a file
void read_field(Field& field, std::string filename,
                const ParallelData parallel)
{
    std::ifstream file;
    file.open(filename);
    // Read the header
    std::string line, comment;
    std::getline(file, line);
    int nx_full, ny_full;
    std::stringstream(line) >> comment >> nx_full >> ny_full;

    field.setup(nx_full, ny_full, parallel);

    auto full = Matrix<double> (nx_full, ny_full);

    // Inner region (no boundaries)
    auto inner = Matrix<double> (field.nx, field.ny);

    if (0 == parallel.rank) {
        // Read the full array
        for (int i = 0; i < nx_full; i++)
            for (int j = 0; j < ny_full; j++)
                file >> full(i, j);

        for (int i = 0; i < field.nx; i++)
            for (int j = 0; j < field.ny; j++)
                inner(i, j) = full(i, j);
    }
    file.close();

    int numelem = field.nx * field.ny;
    MPI_Scatter(full.data(), numelem, MPI_DOUBLE, inner.data(), numelem, MPI_DOUBLE,0, MPI_COMM_WORLD);


    // Copy to the array containing also boundaries
    for (int i = 0; i < field.nx; i++)
        for (int j = 0; j < field.ny; j++)
             field(i + 1, j + 1) = inner(i, j);

    // Set the boundary values
    for (int i = 0; i < field.nx + 2; i++) {
        // left boundary
        field(i, 0) = field(i, 1);
        // right boundary
        field(i, field.ny + 1) = field(i, field.ny);
    }
    for (int j = 0; j < field.ny + 2; j++) {
        // top boundary
        field.temperature(0, j) = field(1, j);
        // bottom boundary
        field.temperature(field.nx + 1, j) = field.temperature(field.nx, j);
    }

}
