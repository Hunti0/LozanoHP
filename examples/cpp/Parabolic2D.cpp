#include <iostream>
#include "mole.h"

int main() {
    int k = 2;
    double a = 0, b = 1, tf = 1;
    int m = 2 * k + 1;
    double dx = (b - a) / m;
    int N = m * m;
    double dt = tf / ceil((4 * tf) / (dx * dx));

    Gradient Gx(k, m, dx);
    Gradient Gy(k, m, dx);
    Divergence Dx(k, m, dx);
    Divergence Dy(k, m, dx);
    identity I(m);

    mat Gx2D = kron(I, Gx);
    mat Gy2D = kron(Gy, I);
    mat Dx2D = kron(I, Dx);
    mat Dy2D = kron(Dy, I);

    // Solution vector
    vec u(N, fill::zeros);

    // Anisotropic diffusion coefficients
    vec alpha_x(N);
    vec alpha_y(N);
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) {
            int idx = i + j * m;
            double x = a + (i + 0.5) * dx;
            double y = a + (j + 0.5) * dx;

            // Example: anisotropic diffusion increasing in x and y
            alpha_x(idx) = 1.0 + x; // e.g., alpha_x ∈ [1,2]
            alpha_y(idx) = 1.0 + y; // e.g., alpha_y ∈ [1,2]
        }
    }

    double t = 0;
    vec k1(N);

    while (t < tf) {
        vec grad_x = Gx2D * u;
        vec grad_y = Gy2D * u;

        // Multiply element-wise with anisotropic diffusion
        vec flux_x = -alpha_x % grad_x;
        vec flux_y = -alpha_y % grad_y;

        vec div_flux_x = Dx2D * flux_x;
        vec div_flux_y = Dy2D * flux_y;

        k1 = div_flux_x + div_flux_y;

        u = u + dt * k1;
        t += dt;
    }

    mat u2D(m, m);
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < m; ++i)
            u2D(i, j) = u(i + j * m);

    std::cout << "Anisotropic solution u(x,y,t):\n";
    std::cout << u2D << std::endl;

    return 0;
}
