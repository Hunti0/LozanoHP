#include <iostream>
#include <vector>
#include <cmath>
#include <memory>  
#include <mole.h>

using namespace std;

// Class for Flux Operators
class FluxOperator {
public:
    virtual ~FluxOperator() = default;  // Virtual destructor ensures proper cleanup of derived classes

    // Pure virtual method to compute flux; must be overridden
    virtual void computeFlux(
        const vector<double>& du_dx,
        const vector<double>& du_dy,
        const vector<double>& du_dz,
        vector<double>& Fx,
        vector<double>& Fy,
        vector<double>& Fz
    ) const = 0;
};

// Concrete implementation of an anisotropic tensor flux
class AnisotropicFlux : public FluxOperator {
private:
    const vector<double>& D11;
    const vector<double>& D12;
    const vector<double>& D13;
    const vector<double>& D21;
    const vector<double>& D22;
    const vector<double>& D23;
    const vector<double>& D31;
    const vector<double>& D32;
    const vector<double>& D33;

public:
    AnisotropicFlux(
        const vector<double>& D11, const vector<double>& D12, const vector<double>& D13,
        const vector<double>& D21, const vector<double>& D22, const vector<double>& D23,
        const vector<double>& D31, const vector<double>& D32, const vector<double>& D33
    ) : D11(D11), D12(D12), D13(D13), D21(D21), D22(D22), D23(D23), D31(D31), D32(D32), D33(D33) {}

    void computeFlux(
        const vector<double>& du_dx,
        const vector<double>& du_dy,
        const vector<double>& du_dz,
        vector<double>& Fx,
        vector<double>& Fy,
        vector<double>& Fz
    ) const override {
        size_t n = Fx.size();
        for (size_t i = 0; i < n; ++i) {
            Fx[i] = -(D11[i] * du_dx[i] + D12[i] * du_dy[i] + D13[i] * du_dz[i]);
            Fy[i] = -(D21[i] * du_dx[i] + D22[i] * du_dy[i] + D23[i] * du_dz[i]);
            Fz[i] = -(D31[i] * du_dx[i] + D32[i] * du_dy[i] + D33[i] * du_dz[i]);
        }
    }
};

int main() {
    // ------------------ Setup Grid ------------------
    int k = 2; // Order of Accuracy
    int nx = 20, ny = 20, nz = 20;  // Number of grid points
    double xmin = 0.0, xmax = 1.0;
    double ymin = 0.0, ymax = 1.0;
    double zmin = 0.0, zmax = 1.0;

    vector<double> x(nx*ny*nz), y(nx*ny*nz), z(nx*ny*nz);

    for (int kk = 0; kk < nz; ++kk) {
        double zeta = static_cast<double>(kk) / (nz - 1);
        for (int j = 0; j < ny; ++j) {
            double eta = static_cast<double>(j) / (ny - 1);
            for (int i = 0; i < nx; ++i) {
                double xi = static_cast<double>(i) / (nx - 1);
                int idx = kk * nx * ny + j * nx + i;

                x[idx] = (1.0 + 0.1 * sin(2 * M_PI * zeta)) * (xmin + xi * (xmax - xmin));
                y[idx] = (1.0 + 0.1 * cos(2 * M_PI * xi)) * (ymin + eta * (ymax - ymin));
                z[idx] = (1.0 + 0.1 * sin(2 * M_PI * eta)) * (zmin + zeta * (zmax - zmin));
            }
        }
    }

    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);
    double dz = (zmax - zmin) / (nz - 1);

    // ------------------ Initial Conditions ------------------
    vector<double> u(nx*ny*nz);
    double x0 = 0.5, y0 = 0.5, z0 = 0.5, sigma = 0.1;

    for (int idx = 0; idx < nx*ny*nz; ++idx) {
        double dx0 = x[idx] - x0;
        double dy0 = y[idx] - y0;
        double dz0 = z[idx] - z0;
        u[idx] = exp(-(dx0*dx0 + dy0*dy0 + dz0*dz0) / (2 * sigma * sigma));
    }

    // ------------------ Setup Flux and Operators ------------------
    vector<double> D11(nx*ny*nz, 1.0);
    vector<double> D12(nx*ny*nz, 0.1);
    vector<double> D13(nx*ny*nz, 0.05);
    vector<double> D21 = D12;
    vector<double> D22(nx*ny*nz, 2.0);
    vector<double> D23(nx*ny*nz, 0.05);
    vector<double> D31 = D13;
    vector<double> D32 = D23;
    vector<double> D33(nx*ny*nz, 0.5);

    unique_ptr<FluxOperator> flux = make_unique<AnisotropicFlux>(
        D11, D12, D13, D21, D22, D23, D31, D32, D33
    );

    Gradient G(k, nx, ny, nz, dx, dy, dz);
    Div3D D(k, nx, ny, nz, dx, dy, dz);

    // ------------------ Time-stepping ------------------
    double dt = 0.0005, t_final = 0.05;
    int n = nx * ny * nz;

    vector<double> Fx(n), Fy(n), Fz(n);

    for (double t = 0; t < t_final; t += dt) {
        // 1. Compute Gradient
        vector<double> grad_u = G * u;

        vector<double> du_dx(grad_u.begin(), grad_u.begin() + n - 1);
        vector<double> du_dy(grad_u.begin() + n - 1, grad_u.begin() + 2 * n);
        vector<double> du_dz(grad_u.begin() + 2 * n, grad_u.begin() + 3 * n);

        // 2. Compute Flux
        flux->computeFlux(du_dx, du_dy, du_dz, Fx, Fy, Fz);

        // 3. Combine Flux for Divergence
        vector<double> F(3 * n);
        copy(Fx.begin(), Fx.end(), F.begin());
        copy(Fy.begin(), Fy.end(), F.begin() + n);
        copy(Fz.begin(), Fz.end(), F.begin() + 2 * n);

        // 4. Compute Divergence of Flux
        vector<double> Lu = D * F;

        // 5. Update u
        for (int idx = 0; idx < n; ++idx) {
            u[idx] += dt * Lu[idx];
        }

        // 6. Apply Dirichlet Boundary Conditions
        for (int kk = 0; kk < nz; ++kk) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (i == 0 || i == nx-1 || j == 0 || j == ny-1 || kk == 0 || kk == nz-1) {
                        int idx = kk * nx * ny + j * nx + i;
                        u[idx] = 0.0;
                    }
                }
            }
        }
    }

    return 0;
}
