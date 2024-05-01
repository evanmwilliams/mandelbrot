#pragma once

#include <vector>
#include <complex>

int mandelbrot(const std::complex<double> &c, int max_iterations);

int map_color(int iterations, int max_iterations);

std::vector<int> generate_mandelbrot_set(
    int width,
    int height,
    double x_min,
    double x_max,
    double y_min,
    double y_max,
    int max_iterations);