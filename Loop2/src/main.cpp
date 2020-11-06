#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
        // anti-dependence in this task
        // sending array's sizes
        MPI_Bcast(&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // sending array from 0 to others
        if (rank != 0)
                arr = (double *) calloc (xSize * ySize, sizeof (double));
        MPI_Bcast (arr, xSize * ySize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // tmp array for calculations
        double* arr_new     = (double *) calloc (xSize * ySize, sizeof (double));

        uint32_t size_part = xSize * ySize / size;
        // copying not changing array's vlues
        if(rank == 0)
        {
                for (uint32_t i = 0; i < xSize * ySize; i++)
                        if ( ((i % xSize) < 3) || ((i / xSize) >= (ySize - 1)) )
                                arr_new[i] = arr[i];
        }
        // calculating
        for (uint32_t i = size_part * rank; i < size_part * (rank + 1); i++)
                if ( ((i % xSize) >= 3) && ((i / xSize) < (ySize - 1)) )
                        arr_new[i] = sin(0.00001*arr[i + xSize - 3]);
        // last process calculates rest of data
        if (rank == size - 1)
        {
                for (uint32_t i = size_part * (rank + 1); i < xSize * ySize; i++)
                        if ( ((i % xSize) >= 3) && ((i / xSize) < (ySize - 1)) )
                                arr_new[i] = sin(0.00001*arr[i + xSize - 3]);
        }
        // gathering data
        MPI_Reduce (arr_new, arr, xSize * ySize, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        // free memory
        if(rank)
                free(arr);
        free(arr_new);

        /*
  if (rank == 0 && size > 0)
  {
    for (uint32_t y = 0; y < ySize - 1; y++)
    {
      for (uint32_t x = 3; x < xSize; x++)
      {
        arr[y*xSize + x] = sin(0.00001*arr[(y + 1)*xSize + x - 3]);
      }
    }
  }
        */

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
