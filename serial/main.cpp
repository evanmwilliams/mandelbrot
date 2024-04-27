#include "common.h"
#include <iostream>

// int main(int argc, char** argv)
// {
//   render();
// }

#include <iostream>
#include <fstream>
#include <complex>

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

int main()
{
  // Image dimensions and parameters
  const int width = 800;
  const int height = 800;
  const double x_min = -2.0, x_max = 1.0;
  const double y_min = -1.5, y_max = 1.5;
  const int max_iterations = 1000;

  // Create a PPM file
  std::ofstream image("mandelbrot.ppm");
  image << "P3\n"
        << width << " " << height << "\n255\n";

  // Generate the Mandelbrot set
  for (int y = 0; y < height; ++y)
  {
    double imag = y_min + (y_max - y_min) * y / (height - 1);
    for (int x = 0; x < width; ++x)
    {
      double real = x_min + (x_max - x_min) * x / (width - 1);
      std::complex<double> c(real, imag);
      int iterations = mandelbrot(c, max_iterations);
      int color = map_color(iterations, max_iterations);

      // Output RGB values (grayscale in this example)
      image << color << " " << color << " " << color << " ";
    }
    image << "\n";
  }

  image.close();
  std::cout << "Mandelbrot image created: mandelbrot.ppm\n";
  return 0;
}
