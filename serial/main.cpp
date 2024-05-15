#include "common.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <chrono>

#define HISTOGRAM 1

std::vector<Color> assign_colors()
{
  std::vector<Color> colors;
  for (int i = 0; i < NUM_ROWS; i++)
  {
    // randomly sample a color from the gradient
    int color_idx = rand() % gradient_samples.size();
    std::cout << "Color index: " << color_idx << std::endl;
    Color color = gradient_samples[color_idx];
    colors.push_back(color);
  }
  // make sure adjacent rows are not same color
  for (int i = 1; i < NUM_ROWS - 1; i++)
  {
    while (colors[i] == colors[i + 1] || colors[i] == colors[i - 1])
    {
      colors[i] = gradient_samples[rand() % gradient_samples.size()];
    }
  }

  return colors;
}

int main()
{
  auto gen_start = std::chrono::high_resolution_clock::now();

  // Image dimensions and parameters

  const double x_min = -2.0, x_max = 1.0;
  const double y_min = -1.5, y_max = 1.5;
  const int max_iterations = 1000;

  std::vector<Color> center_colors = assign_colors();
  // std::cout << "Center colors size: " << center_colors.size() << " and colors[0] is " << center_colors[0].r <<  std::endl;
  //  Generate the Mandelbrot set
  std::vector<Color> mandelbrot_set;
  if (HISTOGRAM)
  {
    std::vector<int> mandelbrot = generate_mandelbrot_set_histogram(x_min, x_max, y_min, y_max);
    mandelbrot_set = color_histogram(mandelbrot, center_colors);
  }
  else
  {
    mandelbrot_set = generate_mandelbrot_set(x_min, x_max, y_min, y_max, center_colors);
  }

  auto gen_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = gen_end - gen_start;

  // Create a PPM file
  std::ofstream image("mandelbrot1.ppm");
  image << "P3\n"
        << WIDTH << " " << HEIGHT << "\n255\n";

  for (int y = 0; y < HEIGHT; ++y)
  {
    for (int x = 0; x < WIDTH; ++x)
    {
      Color c = mandelbrot_set[y * WIDTH + x];
      //if (x % 100 == 0 && y % 100 == 0)
      //  std::cout << "Color: " << c.r << " " << c.g << " " << c.b << std::endl;
      image << c.r << " " << c.g << " " << c.b << " ";
      // image << 255 << " " << 0 << " " << 0 << " ";
    }
    image << "\n";
  }

  image.close();

  std::cout << "Mandelbrot image created: mandelbrot1.ppm\n";
  std::cout << "Time taken to generate image: " << elapsed.count() << " seconds\n";

  return 0;
}
