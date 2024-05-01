#pragma once

#include <iostream>
#include <vector>
#include <complex>

#define MAX_ITER 1000
#define WIDTH 1200
#define HEIGHT 1200
#define NUM_ROWS 4

struct Color {
    int r;
    int g;
    int b;

    bool operator==(const Color &other) const {
        return r == other.r && g == other.g && b == other.b;
    }
};

int mandelbrot(const std::complex<double> &c, int max_iterations);

int map_color(const std::vector<Color>& center_colors, int iterations, int max_iterations);

std::vector<Color> generate_mandelbrot_set(
    double x_min,
    double x_max,
    double y_min,
    double y_max,
    int max_iterations,
    const std::vector<Color>& center_colors);


/*
<palette>
  <color name="Cerulean" hex="0081a7" r="0" g="129" b="167" />
  <color name="Turquoise" hex="45f0df" r="69" g="240" b="223" />
  <color name="Mountbatten pink" hex="9e7b9b" r="158" g="123" b="155" />
  <color name="Wisteria" hex="cb9cf2" r="203" g="156" b="242" />
  <color name="Space cadet" hex="111d4a" r="17" g="29" b="74" />
</palette>
*/
static std::vector<Color> gradient_samples = {
    {0, 129, 167},
    {69, 240, 223},
    {158, 123, 155},
    {203, 156, 242},
    {17, 29, 74}
};
