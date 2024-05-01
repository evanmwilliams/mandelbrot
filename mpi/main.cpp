#include "common.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <chrono>

int main()
{
  auto gen_start = std::chrono::high_resolution_clock::now();

  // Image dimensions and parameters
  const int width = 1200;
  const int height = 1200;
  const double x_min = -2.0, x_max = 1.0;
  const double y_min = -1.5, y_max = 1.5;
  const int max_iterations = 1000;

     // Init MPI
    int num_procs, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // Calculate the number of rows each process will calculate
  int rows_per_process = height / num_procs;
  int start_row = rank * rows_per_process;
  int end_row = (rank == num_procs - 1) ? height : start_row + rows_per_process;

  int proc_y_min = y_min + (y_max - y_min) * start_row / height;
  int proc_y_max = y_min + (y_max - y_min) * end_row / height;

  int proc_height = end_row - start_row;

  // Generate the Mandelbrot set
  std::vector<int> mandelbrot_set = generate_mandelbrot_set(width, proc_height, x_min, x_max, y_min, y_max, max_iterations);

  auto gen_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = gen_end - gen_start;

  std::cout << "Time taken to generate mandelbrot set: " << elapsed.count() << " seconds\n";
  MPI_Finalize();

  // Create a PPM file
  std::ofstream image("mandelbrot_mpi.ppm");
  image << "P3\n"
        << width << " " << height << "\n255\n";

  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      image << mandelbrot_set[y * width + x] << " " << mandelbrot_set[y * width + x] << " " << mandelbrot_set[y * width + x] << " ";
    }
    image << "\n";
  }

  image.close();

  std::cout << "Mandelbrot image created: mandelbrot.ppm\n";
  std::cout << "Time taken to generate image: " << elapsed.count() << " seconds\n";

  return 0;
}
