#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

double acceleration(double t)
{
  return sin(t);
}

void calc(double* trace, uint32_t traceSize, double t0, double dt, double y0, double y1, int rank, int size)
{
        MPI_Status status;
        MPI_Bcast (&traceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast (&t0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast (&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast (&y0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast (&y1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        uint32_t first     = 0;
        uint32_t last      = 0;
        uint32_t part_size = 0;

        if (rank != size - 1)
        {
                first     = traceSize * rank / size;
                last      = traceSize * (rank + 1) / size;
                part_size = last - first;
        }
        else
        {
                first     = traceSize * rank / size;
                last      = traceSize;
                part_size = last - first;
        }

        double v0 = 0.0;
        double start_y = y0;
        double start_v = 0.0;

        double* part_trace = (double* ) malloc(part_size * sizeof(double));

        part_trace[0] = start_y;
        part_trace[1] = start_y + dt * start_v;
        t0 = t0 + (rank * dt * traceSize)/ size;
        for (uint32_t i = 2; i < part_size; i++)
        {
                part_trace[i] = dt*dt*acceleration(t0 + (i - 1)*dt) + 2*part_trace[i - 1] - part_trace[i - 2];
        }
        double end_y = part_trace[part_size - 1];
        double end_v = (part_trace[part_size - 1] - part_trace[part_size - 2])/dt;

        double prev_y = 0.0;
        double prev_v = 0.0;

        if (size == 1)
        {
                v0      = (y1 - end_y) / (dt * traceSize);
                start_v = v0;
        }
        else
        {
                if (rank > 0)
                {
                        MPI_Recv (&prev_y, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
                        MPI_Recv (&prev_v, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
                        start_y = prev_y;
                        start_v = prev_v;
                        end_y += prev_y + prev_v * dt * traceSize/size;
                        end_v += prev_v;
                }
                if (rank != size - 1)
                {
                        MPI_Send (&end_y, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                        MPI_Send (&end_v, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
                }
                if (rank == size - 1)
                {
                        MPI_Send (&end_y, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                }
                if (rank == 0)
                {
                        MPI_Recv (&prev_y, 1, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD, &status);
                        v0 = (y1 - prev_y) / (dt * traceSize);
                }
                MPI_Bcast (&v0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                start_y += (v0 * dt * traceSize * rank)/size;
                start_v += v0;

        }


        part_trace[0] = start_y;
        part_trace[1] = start_y + dt * start_v;
        for (uint32_t i = 2; i < part_size; i++)
        {
                part_trace[i] = dt * dt * acceleration (t0 + (i - 1) * dt)  + 2 * part_trace[i - 1] - part_trace[i - 2];
        }
        if (rank == 0)
        {
                memcpy (trace, part_trace, part_size * sizeof (double));
                for (int i = 1; i < size; i++)
                {
                        uint32_t first = traceSize * i / size;
                        uint32_t last = traceSize * (i + 1) / size;
                        uint32_t len = last - first;
                        MPI_Recv (trace + first, len, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                }
        }
        if (rank > 0)
                MPI_Send (part_trace, part_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

        free (part_trace);

}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, status = 0;
  uint32_t traceSize = 0;
  double t0 = 0, t1 = 0, dt = 0, y0 = 0, y1 = 0;
  double* trace = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    // Check arguments
    if (argc != 3)
    {
      std::cout << "[Error] Usage <inputfile> <output file>\n";
      status = 1;
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Prepare input file
    std::ifstream input(argv[1]);
    if (!input.is_open())
    {
      std::cout << "[Error] Can't open " << argv[1] << " for write\n";
      status = 1;
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Read arguments from input
    input >> t0 >> t1 >> dt >> y0 >> y1;
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    traceSize = (t1 - t0)/dt;
    trace = new double[traceSize];

    input.close();
  } else {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status != 0)
    {
      return 1;
    }
  }

  calc(trace, traceSize, t0, dt, y0, y1, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete trace;
      return 1;
    }

    for (uint32_t i = 0; i < traceSize; i++)
    {
      output << " " << trace[i];
    }
    output << std::endl;
    output.close();
    delete trace;
  }

  MPI_Finalize();
  return 0;
}
