#include "common.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <chrono>
#include <cstdlib>
#include <cassert>
#include <upcxx/upcxx.hpp>

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

std::vector<Color> generate_mandelbrot_set(double x_min, double x_max, double y_min, double y_max, const std::vector<Color>& center_colors)
{
    std::vector<Color> mandelbrot_set(WIDTH * HEIGHT);
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

Color map_hist_color(float hue, std::vector<Color> palette){
    Color color;
    if (hue < 0.2) {
        color = palette[0];
    } else if (hue < 0.4) {
        color = palette[1];
    } else if (hue < 0.6) {
        color = palette[2];
    } else if (hue < 0.96) {
        color = palette[3];
    } else {
        color = palette[4];
    }
    return color;
}

std::vector<Color> generate_mandelbrot_set_histogram(double x_min, double x_max, double y_min, double y_max, const std::vector<Color>& palette)
{
    int y_start = upcxx::rank_me() * (HEIGHT / upcxx::rank_n());  //upcxx::rank_me() * ((HEIGHT / upcxx::rank_n()) + HEIGHT % upcxx::rank_n());
    int y_end = std::min((upcxx::rank_me() + 1) * (HEIGHT / upcxx::rank_n()), HEIGHT); //std::min((upcxx::rank_me() + 1) * ((HEIGHT / upcxx::rank_n()) + HEIGHT % upcxx::rank_n()), HEIGHT);
    int MY_PROC_HEIGHT = y_end - y_start;
    std::cout << "proc " << upcxx::rank_me() << " y_start " << y_start << " y_end " << y_end << std::endl;

    std::vector<int> mandelbrot_set(WIDTH * MY_PROC_HEIGHT);
    for (int y = y_start; y < y_end; ++y)
    {
        double imag = y_min + (y_max - y_min) * y / (HEIGHT - 1);
        for (int x = 0; x < WIDTH; ++x)
        {
            double real = x_min + (x_max - x_min) * x / (WIDTH - 1);
            std::complex<double> c(real, imag);
            int iterations = mandelbrot(c);
            mandelbrot_set[(y - y_start) * WIDTH + x] = iterations;
        }
    }

    std::vector<int> iterations_pp(MAX_ITER + 1, 0);
    float total_before_bail = 0.0;
    for (int i = 0; i < WIDTH * MY_PROC_HEIGHT; ++i) {
        iterations_pp[mandelbrot_set[i]]++;
        if(mandelbrot_set[i] < MAX_ITER) {
            total_before_bail+=1.0;
        }
    }
    upcxx::barrier();

    upcxx::dist_object<std::vector<int>> final_iterations_pp({});
    upcxx::dist_object<float> final_total_before_bail(total_before_bail);
    for(int i = 0; i < iterations_pp.size(); i++) {
        final_iterations_pp->push_back(iterations_pp[i]);
    }
    *final_total_before_bail = total_before_bail;
    std::vector<upcxx::future<>> futs;
    for(int i = 0; i < upcxx::rank_n(); i++) {
        if(i != upcxx::rank_me()) {
            futs.push_back(upcxx::rpc(i, 
                [](upcxx::dist_object<std::vector<int>>& final_iterations_pp, 
                    const std::vector<int>& iterations_pp,  
                    upcxx::dist_object<float>& final_total_before_bail,                 
                    const float total) {
                    for(size_t i = 0; i < iterations_pp.size(); i++) {
                        (*final_iterations_pp)[i] += iterations_pp[i];
                    }
                    *final_total_before_bail += total;
                }, final_iterations_pp, iterations_pp, final_total_before_bail, total_before_bail)
            );
        }
    }
    for(auto& fut : futs) {
        fut.wait();
    }

    std::vector<Color> mandelbrot_colored(WIDTH * MY_PROC_HEIGHT);

    std::vector<int> local_iterations_pp = *final_iterations_pp;
    total_before_bail = *final_total_before_bail;
    for (int y = 0; y < MY_PROC_HEIGHT ; y++) {
        for (int x = 0; x < WIDTH; x++){
            int iteration = mandelbrot_set[y * WIDTH + x];
            float hue_val = 0.0;
            for (int i = 0; i <= iteration; i++) {
                hue_val += (float)((local_iterations_pp)[i]) / total_before_bail;
            }
            mandelbrot_colored[y * WIDTH + x] = map_hist_color(hue_val, palette);
        }
    }

    // send all to rank 0
    upcxx::dist_object<std::vector<Color>> mandelbrot_colored_all(std::vector<Color>(WIDTH * HEIGHT));
    std::vector<upcxx::future<>> futs2;

    for(int y = y_start; y < y_end; y++) {
        for(int x = 0; x < WIDTH; x++) {
            (*mandelbrot_colored_all)[y * WIDTH + x] = mandelbrot_colored[(y - y_start) * WIDTH + x];
        }
    }

    for(int i = 0; i < upcxx::rank_n(); i++) {
        if(i != upcxx::rank_me()) {
            futs2.push_back(upcxx::rpc(i, 
                        [](upcxx::dist_object<std::vector<Color>>& mandelbrot_colored_all, const std::vector<Color>& mandelbrot_colored, int y_start, int y_end) {
                            for(int y = y_start; y < y_end; y++) {
                                for(int x = 0; x < WIDTH; x++) {
                                    (*mandelbrot_colored_all)[y * WIDTH + x] = mandelbrot_colored[(y - y_start) * WIDTH + x];
                                }
                            }
                        }, mandelbrot_colored_all, mandelbrot_colored, y_start, y_end));
        }
    }
    for(auto& fut : futs2) {
        fut.wait();
    }
    upcxx::barrier();

    return *mandelbrot_colored_all;

}
