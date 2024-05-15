#include "common.h"
#include <vector>
#include <complex>
#include <omp.h> // OpenMP header for parallelism

int map_pixel_to_row(int pixel_y) {
    return pixel_y * NUM_ROWS / HEIGHT;
}

// Function to determine if a point is in the Mandelbrot set
int mandelbrot(const std::complex<double> &c)
{
  std::complex<double> z = 0;
  int iterations = 0;
  while (std::norm(z) <= 4 && iterations < MAX_ITER)
  {
    z = z * z + c;
    ++iterations;
  }
  return iterations;
}

int height_of_proc(int proc) {
    int height_base = HEIGHT / NUM_ROWS;
    return (proc == NUM_ROWS - 1) ? HEIGHT - height_base * (NUM_ROWS - 1) : height_base;
}

Color lerp(Color top, Color bot, int t, int max_t) {
    return {
        top.r + (bot.r - top.r) * t / max_t,
        top.g + (bot.g - top.g) * t / max_t,
        top.b + (bot.b - top.b) * t / max_t
    };
}

// Function to map iterations to grayscale color
Color map_color(const std::vector<Color>& center_colors, int iterations, int y_pixel)
{
    if(iterations < MAX_ITER) {
        return {0, 0, 0};
    }
    //std::cout << "Mapping color for pixel " << y_pixel << std::endl;
    int proc_row = map_pixel_to_row(y_pixel);
    int above_row = (proc_row == 0) ? 0 : proc_row - 1;
    int below_row = (proc_row == NUM_ROWS - 1) ? NUM_ROWS - 1 : proc_row + 1;

    int start_pixel_y = proc_row * HEIGHT / NUM_ROWS;
    int center_pixel_y = start_pixel_y + height_of_proc(proc_row) / 2;
    int dist_between_proc = height_of_proc(proc_row);

    //std::cout << "Center color size: " << center_colors.size() << " indexing with i " << proc_row << std::endl;

    if (y_pixel < center_pixel_y) {
        //std::cout << "Pixel " << y_pixel << " is above center " << center_pixel_y << std::endl;
        int dist_to_center = center_pixel_y - y_pixel;
        int dist_to_above = dist_between_proc - dist_to_center;
        if(above_row == proc_row) {
            // fade to white
            return lerp({255, 255, 255}, center_colors[proc_row], dist_to_above, dist_between_proc);
        } else {
            // fade to above proc
            return lerp(center_colors[above_row], center_colors[proc_row], dist_to_above, dist_between_proc);
        }
    } else if (y_pixel == center_pixel_y) {
        return center_colors[proc_row];
    } else {
        int dist_to_center = y_pixel - center_pixel_y;
        int dist_to_above = dist_to_center;
        //std::cout << "Dist to center: " << dist_to_center << std::endl;
        if(below_row == proc_row) {
            // fade to white
            return lerp(center_colors[proc_row], {255, 255, 255}, dist_to_center, dist_between_proc);
        } else {
            // fade to below proc
            return lerp(center_colors[proc_row], center_colors[below_row], dist_to_center, dist_between_proc);
        }
    }
}

// Parallelized function to generate the Mandelbrot set
std::vector<Color> generate_mandelbrot_set(
      double x_min, 
      double x_max, 
      double y_min, 
      double y_max, 
      const std::vector<Color>& center_colors)
{
  std::vector<Color> mandelbrot_set(WIDTH * HEIGHT);

#pragma omp parallel for
  for (int y = 0; y < HEIGHT; ++y)
  {
    double imag = y_min + (y_max - y_min) * y / (HEIGHT - 1);
    for (int x = 0; x < WIDTH; ++x)
    {
      double real = x_min + (x_max - x_min) * x / (WIDTH - 1);
      std::complex<double> c(real, imag);
      int iterations = mandelbrot(c);
      mandelbrot_set[y * WIDTH + x] = map_color(center_colors, iterations, y);
    }
  }
  return mandelbrot_set;
}

std::vector<int> generate_mandelbrot_set_histogram(double x_min, double x_max, double y_min, double y_max)
{
    std::vector<int> mandelbrot_set(WIDTH * HEIGHT);
    for (int y = 0; y < HEIGHT; ++y)
    {
        double imag = y_min + (y_max - y_min) * y / (HEIGHT - 1);
        for (int x = 0; x < WIDTH; ++x)
        {
            double real = x_min + (x_max - x_min) * x / (WIDTH - 1);
            std::complex<double> c(real, imag);
            int iterations = mandelbrot(c);
            mandelbrot_set[y * WIDTH + x] = iterations;
        }
    }
    return mandelbrot_set;
}

std::vector<Color> color_histogram(const std::vector<int>& mandelbrot_set, const std::vector<Color>& palette) {
    std::vector<int> iterations_pp(MAX_ITER + 1, 0);
    float total_before_bail = 0.0;
    for (int i = 0; i < WIDTH * HEIGHT; ++i) {
        iterations_pp[mandelbrot_set[i]]++;
        if(mandelbrot_set[i] < MAX_ITER) {
            total_before_bail+=1.0;
        }
    }
    std::vector<std::vector<float>> hue(WIDTH, std::vector<float>(HEIGHT, 0.0));
    for (int x = 0; x < WIDTH; x++){
        for (int y = 0; y < HEIGHT; y++) {
            int iteration = mandelbrot_set[y * WIDTH + x];
            for (int i = 0; i <= iteration; i++) {
                hue[x][y] += iterations_pp[i] / total_before_bail;
            }
        }
    }

    std::vector<Color> mandelbrot_colored(WIDTH * HEIGHT);
    for (int m = 0; m < WIDTH; m++) {
        for (int n = 0; n < HEIGHT; n++) {
            Color color;
            if (hue[m][n] < 0.2) {
                color = palette[0];
            } else if (hue[m][n] < 0.4) {
                color = palette[1];
            } else if (hue[m][n] < 0.6) {
                color = palette[2];
            } else if (hue[m][n] < 0.96) {
                color = palette[3];
            } else {
                color = palette[4];
            }
            mandelbrot_colored[n * WIDTH + m] = color;
        }
    }
    return mandelbrot_colored;
}

