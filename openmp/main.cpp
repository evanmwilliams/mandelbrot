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
  std::vector<Color> mandelbrot_set = generate_mandelbrot_set(x_min, x_max, y_min, y_max, gradient_samples);

  auto gen_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = gen_end - gen_start;

  std::cout << "Time taken to generate mandelbrot set: " << elapsed.count() << " seconds\n";

  // Create a PPM file
  std::ofstream image("mandelbrot1.ppm");
  image << "P3\n"
        << WIDTH << " " << HEIGHT << "\n255\n";

  for (int y = 0; y < HEIGHT; ++y)
  {
    for (int x = 0; x < WIDTH; ++x)
    {
      Color c = mandelbrot_set[y * WIDTH + x];
      image << c.r << " " << c.g << " " << c.b << " ";
    }
    image << "\n";
  }

  image.close();

  std::cout << "Mandelbrot image created: mandelbrot.ppm\n";
  std::cout << "Time taken to generate image: " << elapsed.count() << " seconds\n";

  return 0;
}
