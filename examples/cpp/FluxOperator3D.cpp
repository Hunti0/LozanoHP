#include <iostream>
#include <vector>
#include <cmath>
#include <memory>  
#include <mole.h>

using namespace std;

//class for Flux Operators
class FluxOperator {
public:
    virtual ~FluxOperator() = default;  // Virtual destructor ensures proper cleanup of derived classes

    // Pure virtual method to compute flux; this will be implemented by derived classes
    virtual void computeFlux(
        const vector<double>& du_dx,    // Partial derivative of u w.r.t. x
        const vector<double>& du_dy,    // Partial derivative of u w.r.t. y
        const vector<double>& du_dz,    // Partial derivative of u w.r.t. z
        vector<double>& Fx,             // Flux in the x direction
        vector<double>& Fy,             // Flux in the y direction
        vector<double>& Fz              // Flux in the z direction
    ) const = 0;  // Pure virtual function; must be overridden in derived classes
};

// Concrete implementation of an anisotropic tensor flux
class AnisotropicFlux : public FluxOperator {
private:
    const vector<double>& D11;  // Diffusivity in the xx direction
    const vector<double>& D12;  // Diffusivity in the xy direction
    const vector<double>& D13;  // Diffusivity in the xz direction
    const vector<double>& D21;  // Diffusivity in the yx direction (equal to D12 in symmetry)
    const vector<double>& D22;  // Diffusivity in the yy direction
    const vector<double>& D23;  // Diffusivity in the yz direction
    const vector<double>& D31;  // Diffusivity in the zx direction (equal to D13 in symmetry)
    const vector<double>& D32;  // Diffusivity in the zy direction (equal to D23 in symmetry)
    const vector<double>& D33;  // Diffusivity in the zz direction

public:
    // Constructor to initialize the diffusivity tensors
    AnisotropicFlux(
        const vector<double>& D11, const vector<double>& D12, const vector<double>& D13,
        const vector<double>& D21, const vector<double>& D22, const vector<double>& D23,
        const vector<double>& D31, const vector<double>& D32, const vector<double>& D33
    ) : D11(D11), D12(D12), D13(D13), D21(D21), D22(D22), D23(D23), D31(D31), D32(D32), D33(D33) {}

    // Override the computeFlux method to implement the anisotropic flux calculation
    void computeFlux(
        const vector<double>& du_dx,  // x-component of the gradient
        const vector<double>& du_dy,  // y-component of the gradient
        const vector<double>& du_dz,  // z-component of the gradient
        vector<double>& Fx,           // Flux in the x direction (output)
        vector<double>& Fy,           // Flux in the y direction (output)
        vector<double>& Fz            // Flux in the z direction (output)
    ) const override {
        size_t n = Fx.size();  // The number of grid points (same for Fx, Fy, Fz)

        // Loop through each grid point and calculate the flux components
        for (size_t i = 0; i < n; ++i) {
            // Flux in x-direction
            Fx[i] = -(D11[i]*du_dx[i] + D12[i]*du_dy[i] + D13[i]*du_dz[i]);

            // Flux in y-direction
            Fy[i] = -(D21[i]*du_dx[i] + D22[i]*du_dy[i] + D23[i]*du_dz[i]);

            // Flux in z-direction
            Fz[i] = -(D31[i]*du_dx[i] + D32[i]*du_dy[i] + D33[i]*du_dz[i]);
        }
    }
};

int main() {

    // Grid parameters
    int nx = 20, ny = 20, nz = 20;
    double xmin = 0.0, xmax = 1.0;
    double ymin = 0.0, ymax = 1.0;
    double zmin = 0.0, zmax = 1.0;

    // Allocate vectors to hold the curvilinear grid points (x, y, z) in 3D space
    vector<double> x(nx*ny*nz), y(nx*ny*nz), z(nx*ny*nz);

    // Loop to generate curvilinear grid points with some sinusoidal variation for complexity
    for (int k = 0; k < nz; ++k) {
        double zeta = static_cast<double>(k) / (nz - 1);  // Z-coordinate (normalized)
        for (int j = 0; j < ny; ++j) {
            double eta = static_cast<double>(j) / (ny - 1);  // Y-coordinate (normalized)
            for (int i = 0; i < nx; ++i) {
                double xi = static_cast<double>(i) / (nx - 1);  // X-coordinate (normalized)
                int idx = k * nx * ny + j * nx + i;  // Index into the flattened grid vector

                // Modify coordinates with sinusoidal variations (for curvilinear effect)
                x[idx] = (1.0 + 0.1 * sin(2 * M_PI * zeta)) * (xmin + xi * (xmax - xmin));
                y[idx] = (1.0 + 0.1 * cos(2 * M_PI * xi)) * (ymin + eta * (ymax - ymin));
                z[idx] = (1.0 + 0.1 * sin(2 * M_PI * eta)) * (zmin + zeta * (zmax - zmin));
            }
        }
    }

    // Define grid spacing in each direction
    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);
    double dz = (zmax - zmin) / (nz - 1);

    // Create gradient and divergence operators 
     Gradient G(k,nx,ny,dx,dy);  // Gradient operator (3D)
    Div3D D(nx, ny, nz, dx, dy, dz);   // Divergence operator (3D)

    // Define the anisotropic diffusion tensor components
    vector<double> D11(nx*ny*nz, 1.0);  // Diffusivity in x-direction (Dxx)
    vector<double> D12(nx*ny*nz, 0.1);  // Diffusivity in xy-direction (Dxy)
    vector<double> D13(nx*ny*nz, 0.05); // Diffusivity in xz-direction (Dxz)
    vector<double> D21 = D12;           // Symmetric component (Dyx = Dxy)
    vector<double> D22(nx*ny*nz, 2.0);  // Diffusivity in y-direction (Dyy)
    vector<double> D23(nx*ny*nz, 0.05); // Diffusivity in yz-direction (Dyz)
    vector<double> D31 = D13;           // Symmetric component (Dzx = Dxz)
    vector<double> D32 = D23;           // Symmetric component (Dzy = Dyz)
    vector<double> D33(nx*ny*nz, 0.5);  // Diffusivity in z-direction (Dzz)

    // Abstract flux object - create an instance of the anisotropic tensor flux
    unique_ptr<FluxComputer> flux = make_unique<AnisotropicTensorFlux>(
        D11, D12, D13, D21, D22, D23, D31, D32, D33
    );

    // Initial condition: Gaussian pulse centered at (0.5, 0.5, 0.5)
    vector<double> u(nx*ny*nz);
    double x0 = 0.5, y0 = 0.5, z0 = 0.5, sigma = 0.1;
    for (int idx = 0; idx < nx*ny*nz; ++idx) {
        double dx = x[idx] - x0, dy = y[idx] - y0, dz = z[idx] - z0;
        u[idx] = exp(-(dx*dx + dy*dy + dz*dz) / (2 * sigma * sigma));
    }

    // Time-stepping parameters
    double dt = 0.0005, t_final = 0.05;
    int n = nx * ny * nz;  // Total number of grid points

    // Arrays for flux components in x, y, and z directions
    vector<double> Fx(n), Fy(n), Fz(n);

    // Time loop for simulation
    for (double t = 0; t < t_final; t += dt) {
        // Compute the gradient of u (∇u)
        vector<double> grad_u = G * u;

        // Extract individual components of the gradient (∂u/∂x, ∂u/∂y, ∂u/∂z)
        vector<double> du_dx(grad_u.begin(), grad_u.begin() + n);
        vector<double> du_dy(grad_u.begin() + n, grad_u.begin() + 2 * n);
        vector<double> du_dz(grad_u.begin() + 2 * n, grad_u.begin() + 3 * n);

        // Call the flux operator to compute the fluxes
        flux->computeFlux(du_dx, du_dy, du_dz, Fx, Fy, Fz);

        // Combine the fluxes into a single flux vector (for divergence calculation)
        vector<double> F(3 * n);
        copy(Fx.begin(), Fx.end(), F.begin());
        copy(Fy.begin(), Fy.end(), F.begin() + n);
        copy(Fz.begin(), Fz.end(), F.begin() + 2 * n);

        // Compute the divergence of the flux (∇·F)
        vector<double> Lu = D * F;

        // Update the solution based on the time step
        for (int i = 0; i < n; ++i) {
            u[i] += dt * Lu[i];
        }

        // Apply Dirichlet boundary conditions (u = 0 at the boundaries)
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1 || k == 0 || k == nz - 1) {
                        int idx = k * nx * ny + j * nx + i;
                        u[idx] = 0.0;  // Set boundary points to 0
                    }
                }
            }
        }
    }

    return 0;  // End of the program
}
