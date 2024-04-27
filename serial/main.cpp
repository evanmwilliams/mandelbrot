#include "common.h"
#include <iostream>

// int main(int argc, char** argv)
// {
//   render();
// }

#include <iostream>
#include <fstream>
#include <complex>
#include <chrono>

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

std::vector<int> generate_mandelbrot_set(int width, int height, double x_min, double x_max, double y_min, double y_max, int max_iterations) {
  std::vector<int> mandelbrot_set(width * height);
  for (int y = 0; y < height; ++y)
  {
    double imag = y_min + (y_max - y_min) * y / (height - 1);
    for (int x = 0; x < width; ++x)
    {
      double real = x_min + (x_max - x_min) * x / (width - 1);
      std::complex<double> c(real, imag);
      int iterations = mandelbrot(c, max_iterations);
      mandelbrot_set[y * width + x] = map_color(iterations, max_iterations);
    }
  }
  return mandelbrot_set;

}

int main()
{
  auto gen_start = std::chrono::high_resolution_clock::now();

  // Image dimensions and parameters
  const int width = 1200;
  const int height = 1200;
  const double x_min = -2.0, x_max = 1.0;
  const double y_min = -1.5, y_max = 1.5;
  const int max_iterations = 1000;

  // Generate the Mandelbrot set
  std::vector<int> mandelbrot_set = generate_mandelbrot_set(width, height, x_min, x_max, y_min, y_max, max_iterations);

  auto gen_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = gen_end - gen_start;

  std::cout << "Time taken to generate mandelbrot set: " << elapsed.count() << " seconds\n";

  // Create a PPM file
  std::ofstream image("mandelbrot1.ppm");
  image << "P3\n"
        << width << " " << height << "\n255\n";

  for (int y = 0; y < height; ++y){
    for (int x = 0; x < width; ++x){
      image << mandelbrot_set[y * width + x] << " " << mandelbrot_set[y * width + x] << " " << mandelbrot_set[y * width + x] << " ";
    }
    image << "\n";
  }

  image.close();


  std::cout << "Mandelbrot image created: mandelbrot.ppm\n";
  std::cout << "Time taken to generate image: " << elapsed.count() << " seconds\n";

  return 0;
}
