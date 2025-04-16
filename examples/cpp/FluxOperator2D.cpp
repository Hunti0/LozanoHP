#include <iostream>
#include <vector>
#include <cmath>
#include <mole/mole.h>

using namespace std;

int main() {

    int nx = 50, ny = 50;
    double xmin = 0.0, xmax = 1.0;
    double ymin = 0.0, ymax = 1.0;

    // Curvilinear grid (polar-like transformation)
    vector<double> x(nx*ny), y(nx*ny);
    for (int j = 0; j < ny; ++j) {
        double eta = static_cast<double>(j)/(ny-1);
        for (int i = 0; i < nx; ++i) {
            double xi = static_cast<double>(i)/(nx-1);
            int idx = j*nx + i;
            x[idx] = (1.0 + 0.2*sin(2*M_PI*eta)) * (xmin + xi*(xmax-xmin));
            y[idx] = (1.0 + 0.3*cos(2*M_PI*xi)) * (ymin + eta*(ymax-ymin));
        }
    }


    Grad G(nx, ny, x, y);  // Gradient
    Div D(G);           // Divergence 

    // Anisotropic diffusion tensor
    vector<double> D11(nx*ny, 1.0);  // Dxx
    vector<double> D12(nx*ny, 0.1);  // Dxy
    vector<double> D21(nx*ny, 0.1);  // Dyx
    vector<double> D22(nx*ny, 2.0);  // Dyy

    // Initial condition (Gaussian pulse)
    vector<double> u(nx*ny);
    double x0 = 0.5, y0 = 0.5, sigma = 0.1;
    for (int idx = 0; idx < nx*ny; ++idx) {
        double dx = x[idx] - x0, dy = y[idx] - y0;
        u[idx] = exp(-(dx*dx + dy*dy)/(2*sigma*sigma));
    }

    // Time stepping
    double dt = 0.001, t_final = 0.1;
    for (double t = 0; t < t_final; t += dt) {
        // Compute gradient 
        vector<double> Grad_u = G * u;

        // Split gradient into components
        int n = nx * ny;
        vector<double> du_dx(grad_u.begin(), grad_u.begin() + n);
        vector<double> du_dy(grad_u.begin() + n, grad_u.end());

        // Compute flux F = -K∇u = [Fx; Fy]
        vector<double> Fx(n), Fy(n);
        for (int i = 0; i < n; ++i) {
            Fx[i] = -(D11[i] * du_dx[i] + D12[i] * du_dy[i]);  // Fx = -(D11 ∂u/∂x + D12 ∂u/∂y)
            Fy[i] = -(D21[i] * du_dx[i] + D22[i] * du_dy[i]);  // Fy = -(D21 ∂u/∂x + D22 ∂u/∂y)
        }

        // Combine into single flux vector [Fx; Fy]
        vector<double> F(2 * n);
        copy(Fx.begin(), Fx.end(), F.begin());
        copy(Fy.begin(), Fy.end(), F.begin() + n);

        // Compute divergence of flux ∇·F
        vector<double> Lu = D * F;

        // Update solution
        for (int i = 0; i < n; ++i) {
            u[i] += dt * Lu[i];
        }

        // Boundary conditions (Dirichlet zero)
        for (int i = 0; i < nx; ++i) {
            u[i] = u[(ny-1)*nx + i] = 0.0;  // Bottom/Top
        }
        for (int j = 0; j < ny; ++j) {
            u[j*nx] = u[j*nx + (nx-1)] = 0.0;  // Left/Right
        }
    }

    return 0;
}
