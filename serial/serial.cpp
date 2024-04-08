#include <iostream>

#define MAX_ITER 1000
#define WIDTH 800
#define HEIGHT 600

// Function to compute if a point is in the Mandelbrot set
int compute_point(double x, double y) {
    double z_re = x, z_im = y;
    int n;
    for (n = 0; n < MAX_ITER; ++n) {
        if (z_re * z_re + z_im * z_im > 4)
            break;
        double new_re = z_re * z_re - z_im * z_im;
        double new_im = 2 * z_re * z_im;
        z_re = new_re + x;
        z_im = new_im + y;
    }
    return n;
}

// Main function to render the Mandelbrot set
void render() {
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            double x_scaled = (x - WIDTH / 2.0) * 4.0 / WIDTH;
            double y_scaled = (y - HEIGHT / 2.0) * 4.0 / HEIGHT;
            int iteration = compute_point(x_scaled, y_scaled);
            if (iteration == MAX_ITER) std::cout << ".";
            else std::cout << " ";
        }
        std::cout << std::endl;
    }
}

