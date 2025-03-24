#include <iostream>
#include "mole.h"

int main() {
    int k = 2;
    double a = 0, b = 1, tf = 1;
    int m = 2 * k + 1;
    double dx = (b - a) / m;
    int N = m * m * m;
    double dt = tf / ceil((6 * tf) / (dx * dx)); // tighter for 3D stability

    // 1D mimetic operators
    Gradient G(k, m, dx);
    Divergence D(k, m, dx);
    identity I(m);

    // 3D Kronecker operators
    mat Ix = I, Iy = I, Iz = I;

    mat Gx3D = kron(kron(Iz, Iy), G);
    mat Gy3D = kron(kron(Iz, G), Ix);
    mat Gz3D = kron(kron(G, Iy), Ix);

    mat Dx3D = kron(kron(Iz, Iy), D);
    mat Dy3D = kron(kron(Iz, D), Ix);
    mat Dz3D = kron(kron(D, Iy), Ix);

    // Solution
    vec u(N, fill::zeros);

    // Anisotropic α in x, y, z directions
    vec alpha_x(N), alpha_y(N), alpha_z(N);
    for (int k_idx = 0; k_idx < m; ++k_idx) {
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < m; ++i) {
                int idx = i + j * m + k_idx * m * m;
                double x = a + (i + 0.5) * dx;
                double y = a + (j + 0.5) * dx;
                double z = a + (k_idx + 0.5) * dx;

                alpha_x(idx) = 1.0 + x;
                alpha_y(idx) = 1.0 + y;
                alpha_z(idx) = 1.0 + z;
            }
        }
    }

    double t = 0;
    vec k1(N);

    while (t < tf) {
        vec grad_x = Gx3D * u;
        vec grad_y = Gy3D * u;
        vec grad_z = Gz3D * u;

        vec flux_x = -alpha_x % grad_x;
        vec flux_y = -alpha_y % grad_y;
        vec flux_z = -alpha_z % grad_z;

        vec div_flux_x = Dx3D * flux_x;
        vec div_flux_y = Dy3D * flux_y;
        vec div_flux_z = Dz3D * flux_z;

        k1 = div_flux_x + div_flux_y + div_flux_z;

        u = u + dt * k1;
        t += dt;
    }

    std::cout << "3D anisotropic solution computed!" << std::endl;
    return 0;
}
