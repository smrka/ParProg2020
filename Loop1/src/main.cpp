#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
        // sending array's size
        MPI_Bcast(&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // dividing array into pieces
        int arr_size      = xSize * ySize;
        int arr_part_size = arr_size/size;
        double*  arr_part = (double *)malloc(arr_part_size * sizeof(double));
        MPI_Scatter(arr, arr_part_size, MPI_DOUBLE, arr_part, arr_part_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // calculating
        for(int i = 0; i < arr_part_size; i++)
                arr_part[i] = sin(0.00001 * arr_part[i]);

        // gathering pieces
        MPI_Gather(arr_part, arr_part_size, MPI_DOUBLE, arr, arr_part_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // calculating the rest of input array
        if (rank == 0)
                for(int i = arr_part_size * size; i < arr_size; i++)
                        arr[i] = sin(0.00001 * arr[i]);

        // free memory
        free(arr_part);
}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, buf = 0;
  uint32_t ySize = 0, xSize = 0;
  double* arr = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    // Check arguments
    if (argc != 3)
    {
      std::cout << "[Error] Usage <inputfile> <output file>\n";
      buf = 1;
      MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Prepare input file
    std::ifstream input(argv[1]);
    if (!input.is_open())
    {
      std::cout << "[Error] Can't open " << argv[1] << " for write\n";
      buf = 1;
      MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Read arguments from input
    input >> ySize >> xSize;
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    arr = new double[ySize * xSize];

    for (uint32_t y = 0; y < ySize; y++)
    {
     for (uint32_t x = 0; x < xSize; x++)
      {
        input >> arr[y*xSize + x];
      }
    }
    input.close();
  } else {
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (buf != 0)
    {
      return 1;
    }
  }

  calc(arr, ySize, xSize, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete arr;
      return 1;
    }
    for (uint32_t y = 0; y < ySize; y++)
    {
      for (uint32_t x = 0; x < xSize; x++)
      {
        output << " " << arr[y*xSize + x];
      }
      output << std::endl;
    }
    output.close();
    delete arr;
  }

  MPI_Finalize();
  return 0;
}
