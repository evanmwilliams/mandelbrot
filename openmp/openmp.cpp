#include "common.h"
#include <vector>
#include <complex>
#include <omp.h> // OpenMP header for parallelism

// Function to determine if a point is in the Mandelbrot set
int mandelbrot(const std::complex<double> &c, int max_iterations)
{
  std::complex<double> z = 0;
  int iterations = 0;
  while (std::norm(z) <= 4 && iterations < max_iterations)
  {
    z = z * z + c;
    ++iterations;
  }
  return iterations;
}

// Function to map iterations to grayscale color
int map_color(int iterations, int max_iterations)
{
  return (iterations * 255) / max_iterations;
}

// Parallelized function to generate the Mandelbrot set
std::vector<int> generate_mandelbrot_set(int width, int height, double x_min, double x_max, double y_min, double y_max, int max_iterations)
{
  std::vector<int> mandelbrot_set(width * height);

#pragma omp parallel for
  for (int y = 0; y < height; ++y)
  {
    double imag = y_min + (y_max - y_min) * y / (height - 1);
    for (int x = 0; x < width; ++x)
    {
      double real = x_min + (x_max - x_min) * x / (width - 1);
      std::complex<double> c(real, imag);
      int iterations = mandelbrot(c, max_iterations);
      int color = map_color(iterations, max_iterations);

      mandelbrot_set[y * width + x] = color;
    }
  }

  return mandelbrot_set;
}
