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

  // Generate the Mandelbrot set
  std::vector<int> mandelbrot_set = generate_mandelbrot_set(width, height, x_min, x_max, y_min, y_max, max_iterations);

  auto gen_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = gen_end - gen_start;

  std::cout << "Time taken to generate mandelbrot set: " << elapsed.count() << " seconds\n";

  // Create a PPM file
  std::ofstream image("mandelbrot1.ppm");
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
