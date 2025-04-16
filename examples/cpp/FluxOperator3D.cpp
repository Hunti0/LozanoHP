#include <iostream>
#include <vector>
#include <cmath>
#include <mole/mole.h>

using namespace std;

int main() {

    // Grid
    int nx = 20, ny = 20, nz = 20;
    double xmin = 0.0, xmax = 1.0;
    double ymin = 0.0, ymax = 1.0;
    double zmin = 0.0, zmax = 1.0;

    // Curvilinear grid 
    vector<double> x(nx*ny*nz), y(nx*ny*nz), z(nx*ny*nz);
    for (int k = 0; k < nz; ++k) {
        double zeta = static_cast<double>(k)/(nz-1);
        for (int j = 0; j < ny; ++j) {
            double eta = static_cast<double>(j)/(ny-1);
            for (int i = 0; i < nx; ++i) {
                double xi = static_cast<double>(i)/(nx-1);
                int idx = k*nx*ny + j*nx + i;
                x[idx] = (1.0 + 0.1*sin(2*M_PI*zeta)) * (xmin + xi*(xmax-xmin));
                y[idx] = (1.0 + 0.1*cos(2*M_PI*xi)) * (ymin + eta*(ymax-ymin));
                z[idx] = (1.0 + 0.1*sin(2*M_PI*eta)) * (zmin + zeta*(zmax-zmin));
            }
        }
    }


    Grad3D G(nx, ny, nz, x, y, z);  // 3D Gradient
    Div3D D(G);                   // 3D Divergence

    // Anisotropic diffusion tensor 
    vector<double> D11(nx*ny*nz, 1.0);  // Dxx
    vector<double> D12(nx*ny*nz, 0.1);  // Dxy
    vector<double> D13(nx*ny*nz, 0.05); // Dxz
    vector<double> D21 = D12;           // Dyx
    vector<double> D22(nx*ny*nz, 2.0);  // Dyy
    vector<double> D23(nx*ny*nz, 0.05); // Dyz
    vector<double> D31 = D13;           // Dzx
    vector<double> D32 = D23;           // Dzy
    vector<double> D33(nx*ny*nz, 0.5);  // Dzz

    // Initial condition (Gaussian pulse)
    vector<double> u(nx*ny*nz);
    double x0 = 0.5, y0 = 0.5, z0 = 0.5, sigma = 0.1;
    for (int idx = 0; idx < nx*ny*nz; ++idx) {
        double dx = x[idx] - x0, dy = y[idx] - y0, dz = z[idx] - z0;
        u[idx] = exp(-(dx*dx + dy*dy + dz*dz)/(2*sigma*sigma));
    }

    // Time stepping
    double dt = 0.0005, t_final = 0.05;
    for (double t = 0; t < t_final; t += dt) {
        // Compute gradient ∇u (returns [∂u/∂x; ∂u/∂y; ∂u/∂z])
        vector<double> grad_u = G * u;

        // Split gradient into components
        int n = nx * ny * nz;
        vector<double> du_dx(grad_u.begin(), grad_u.begin() + n);
        vector<double> du_dy(grad_u.begin() + n, grad_u.begin() + 2*n);
        vector<double> du_dz(grad_u.begin() + 2*n, grad_u.end());

        // Compute flux F = -K∇u = [Fx; Fy; Fz]
        vector<double> Fx(n), Fy(n), Fz(n);
        for (int i = 0; i < n; ++i) {
            Fx[i] = -(D11[i] * du_dx[i] + D12[i] * du_dy[i] + D13[i] * du_dz[i]);
            Fy[i] = -(D21[i] * du_dx[i] + D22[i] * du_dy[i] + D23[i] * du_dz[i]);
            Fz[i] = -(D31[i] * du_dx[i] + D32[i] * du_dy[i] + D33[i] * du_dz[i]);
        }

        // Combine into single flux vector [Fx; Fy; Fz]
        vector<double> F(3 * n);
        copy(Fx.begin(), Fx.end(), F.begin());
        copy(Fy.begin(), Fy.end(), F.begin() + n);
        copy(Fz.begin(), Fz.end(), F.begin() + 2*n);

        // Compute divergence of flux ∇·F
        vector<double> Lu = D * F;

        // Update solution
        for (int i = 0; i < n; ++i) {
            u[i] += dt * Lu[i];
        }

        // Boundary conditions (Dirichlet zero)
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (i == 0 || i == nx-1 || j == 0 || j == ny-1 || k == 0 || k == nz-1) {
                        int idx = k*nx*ny + j*nx + i;
                        u[idx] = 0.0;
                    }
                }
            }
        }
    }

    return 0;
}
