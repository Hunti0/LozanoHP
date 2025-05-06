#include <iostream>
#include <vector>
#include <cmath>
#include <memory>  
#include <mole.h>

using namespace std;

// ---------------- Grid Classes ----------------

// Abstract base class
class Grid3D {
public:
    virtual ~Grid3D() = default;
    virtual void generateGrid(vector<double>& x, vector<double>& y, vector<double>& z) const = 0;
};

// Cartesian grid (non-curvilinear)
class CartesianGrid3D : public Grid3D {
private:
    int nx, ny, nz;
    double xmin, xmax, ymin, ymax, zmin, zmax;

public:
    CartesianGrid3D(int nx_, int ny_, int nz_,
                    double xmin_, double xmax_,
                    double ymin_, double ymax_,
                    double zmin_, double zmax_)
        : nx(nx_), ny(ny_), nz(nz_),
          xmin(xmin_), xmax(xmax_),
          ymin(ymin_), ymax(ymax_),
          zmin(zmin_), zmax(zmax_) {}

    void generateGrid(vector<double>& x, vector<double>& y, vector<double>& z) const override {
        x.resize(nx * ny * nz);
        y.resize(nx * ny * nz);
        z.resize(nx * ny * nz);

        for (int i = 0; i < nz; ++i) {
            double zeta = static_cast<double>(i) / (nz - 1);
            for (int j = 0; j < ny; ++j) {
                double eta = static_cast<double>(j) / (ny - 1);
                for (int l = 0; l < nx; ++l) {
                    double xi = static_cast<double>(l) / (nx - 1);
                    int idx = i * nx * ny + j * nx + l;

                    x[idx] = xmin + xi * (xmax - xmin);
                    y[idx] = ymin + eta * (ymax - ymin);
                    z[idx] = zmin + zeta * (zmax - zmin);
                }
            }
        }
    }
};

// Curvilinear grid
class CurvilinearGrid3D : public Grid3D {
private:
    int nx, ny, nz;
    double xmin, xmax, ymin, ymax, zmin, zmax;

public:
    CurvilinearGrid3D(int nx_, int ny_, int nz_,
                      double xmin_, double xmax_,
                      double ymin_, double ymax_,
                      double zmin_, double zmax_)
        : nx(nx_), ny(ny_), nz(nz_),
          xmin(xmin_), xmax(xmax_),
          ymin(ymin_), ymax(ymax_),
          zmin(zmin_), zmax(zmax_) {}

    void generateGrid(vector<double>& x, vector<double>& y, vector<double>& z) const override {
        x.resize(nx * ny * nz);
        y.resize(nx * ny * nz);
        z.resize(nx * ny * nz);

        for (int i = 0; i < nz; ++i) {
            double zeta = static_cast<double>(i) / (nz - 1);
            for (int j = 0; j < ny; ++j) {
                double eta = static_cast<double>(j) / (ny - 1);
                for (int l = 0; l < nx; ++l) {
                    double xi = static_cast<double>(l) / (nx - 1);
                    int idx = i * nx * ny + j * nx + l;

                    x[idx] = (1.0 + 0.1 * sin(2 * M_PI * zeta)) * (xmin + xi * (xmax - xmin));
                    y[idx] = (1.0 + 0.1 * cos(2 * M_PI * xi)) * (ymin + eta * (ymax - ymin));
                    z[idx] = (1.0 + 0.1 * sin(2 * M_PI * eta)) * (zmin + zeta * (zmax - zmin));
                }
            }
        }
    }
};

// ---------------- Flux Operator Classes ----------------

class FluxOperator {
public:
    virtual ~FluxOperator() = default;

    virtual void computeFlux(
        const vector<double>& du_dx,
        const vector<double>& du_dy,
        const vector<double>& du_dz,
        vector<double>& Fx,
        vector<double>& Fy,
        vector<double>& Fz
    ) const = 0;
};

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

// ---------------- Main Program ----------------

int main() {
    // Parameters
    int k = 2;  // Order of Accuracy
    int nx = 20, ny = 20, nz = 20;
    double xmin = 0.0, xmax = 1.0;
    double ymin = 0.0, ymax = 1.0;
    double zmin = 0.0, zmax = 1.0;

    vector<double> x, y, z;

    // Choose grid type
    unique_ptr<Grid3D> grid;
    bool use_curvilinear = false;  // <-- easily switch here

    if (use_curvilinear) {
        grid = make_unique<CurvilinearGrid3D>(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax);
    } else {
        grid = make_unique<CartesianGrid3D>(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax);
    }

    grid->generateGrid(x, y, z);

    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);
    double dz = (zmax - zmin) / (nz - 1);

    int n = nx * ny * nz;

    // Initial condition
    vector<double> u(n);
    double x0 = 0.5, y0 = 0.5, z0 = 0.5, sigma = 0.1;
    for (int idx = 0; idx < n; ++idx) {
        double dx0 = x[idx] - x0;
        double dy0 = y[idx] - y0;
        double dz0 = z[idx] - z0;
        u[idx] = exp(-(dx0*dx0 + dy0*dy0 + dz0*dz0) / (2 * sigma * sigma));
    }

    // Anisotropic diffusion tensor
    vector<double> D11(n, 1.0);
    vector<double> D12(n, 0.1);
    vector<double> D13(n, 0.05);
    vector<double> D21 = D12;
    vector<double> D22(n, 2.0);
    vector<double> D23(n, 0.05);
    vector<double> D31 = D13;
    vector<double> D32 = D23;
    vector<double> D33(n, 0.5);

    unique_ptr<FluxOperator> flux = make_unique<AnisotropicFlux>(
        D11, D12, D13, D21, D22, D23, D31, D32, D33
    );

    Gradient G(k, nx, ny, nz, dx, dy, dz);
    Div3D D(k, nx, ny, nz, dx, dy, dz);

    vector<double> Fx(n), Fy(n), Fz(n);

    // Time-stepping
    double dt = 0.0005, t_final = 0.05;
    for (double t = 0; t < t_final; t += dt) {
        // Gradient
        vector<double> grad_u = G * u;
        vector<double> du_dx(grad_u.begin(), grad_u.begin() + n - 1);
        vector<double> du_dy(grad_u.begin() + n - 1, grad_u.begin() + 2 * n);
        vector<double> du_dz(grad_u.begin() + 2 * n, grad_u.end());

        // Compute Flux
        flux->computeFlux(du_dx, du_dy, du_dz, Fx, Fy, Fz);

        // Combine flux into one vector
        vector<double> F(3 * n);
        copy(Fx.begin(), Fx.end(), F.begin());
        copy(Fy.begin(), Fy.end(), F.begin() + n);
        copy(Fz.begin(), Fz.end(), F.begin() + 2 * n);

        // Divergence of flux
        vector<double> Lu = D * F;

        // Update solution
        for (int i = 0; i < n; ++i) {
            u[i] += dt * Lu[i];
        }

        // Dirichlet boundary conditions
        for (int ii = 0; ii < nz; ++ii) {
            for (int jj = 0; jj < ny; ++jj) {
                for (int ll = 0; ll < nx; ++ll) {
                    if (ll == 0 || ll == nx-1 || jj == 0 || jj == ny-1 || ii == 0 || ii == nz-1) {
                        int idx = ii * nx * ny + jj * nx + ll;
                        u[idx] = 0.0;
                    }
                }
            }
        }
    }

    return 0;
}
