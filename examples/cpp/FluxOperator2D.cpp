#include <iostream>
#include <vector>
#include <cmath>
#include <mole/mole.h>

using namespace std;

int main() {
    
    // Grid parameters
    int nx = 50;
    int ny = 50;
    double xmin = 0.0, xmax = 1.0;
    double ymin = 0.0, ymax = 1.0;
    
    // Create curvilinear grid
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
    
    Grad G(nx, ny, x, y);
    Div D(G);
    
    // Anisotropic diffusion tensor components
    vector<double> D11(nx*ny, 1.0);  // Dxx
    vector<double> D12(nx*ny, 0.1);  // Dxy
    vector<double> D21(nx*ny, 0.1);  // Dyx
    vector<double> D22(nx*ny, 2.0);  // Dyy
    
    // Function to apply the flux operator explicitly: F = -K*∇u
    auto applyFlux = [&](const vector<double>& u) {
        
        vector<double> grad_u = Grad*u;
        
        // The gradient is stored as [du/dx; du/dy]
        int n = nx*ny;
        vector<double> du_dx(grad_u.begin(), grad_u.begin() + n);
        vector<double> du_dy(grad_u.begin() + n, grad_u.end());
        
        // Compute flux components: Fx = -(D11*du_dx + D12*du_dy)
        //                         Fy = -(D21*du_dx + D22*du_dy)
        vector<double> Fx(n), Fy(n);
        for (int i = 0; i < n; ++i) {
            Fx[i] = -(D11[i]*du_dx[i] + D12[i]*du_dy[i]);
            Fy[i] = -(D21[i]*du_dx[i] + D22[i]*du_dy[i]);
        }
        
        // Combine flux components into single vector [Fx; Fy]
        vector<double> F(2*n);
        copy(Fx.begin(), Fx.end(), F.begin());
        copy(Fy.begin(), Fy.end(), F.begin() + n);
        
        return F;
    };
    
    // Initial condition
    vector<double> u(nx*ny);
    double x0 = 0.5, y0 = 0.5;
    double sigma = 0.1;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int idx = j*nx + i;
            double dx = x[idx] - x0;
            double dy = y[idx] - y0;
            u[idx] = exp(-(dx*dx + dy*dy)/(2*sigma*sigma));
        }
    }
    
    // Time stepping
    double dt = 0.001;
    double t_final = 0.1;
    int nsteps = static_cast<int>(t_final/dt);
    
    for (int step = 0; step < nsteps; ++step) {
        // Compute flux F = -K*∇u
        vector<double> F = applyFlux(u);
        
        // Compute divergence of flux ∇·F
        vector<double> Lu = D*F;
        
        // Update solution
        for (int i = 0; i < nx*ny; ++i) {
            u[i] += dt * Lu[i];
        }
        
        // Boundary conditions
        for (int i = 0; i < nx; ++i) {
            u[i] = 0.0;                   // Bottom
            u[(ny-1)*nx + i] = 0.0;       // Top
        }
        for (int j = 0; j < ny; ++j) {
            u[j*nx] = 0.0;                // Left
            u[j*nx + (nx-1)] = 0.0;       // Right
        }
    }
    
    return 0;
}
