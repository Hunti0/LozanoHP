#include <iostream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include <cmath>

#include "mole.h"
#include "Flux3D.h"

namespace {
constexpr double Lx = 1.0;
constexpr double Ly = 1.0;
constexpr double kLambda = 1.0;

inline double U_exact(double x, double y) {
  return std::exp(kLambda * (x + y) / 2.0) / (std::exp(kLambda) - 1.0);
}
inline double Ux_exact(double x, double y) { return 0.5 * kLambda * U_exact(x, y); }
inline double Uy_exact(double x, double y) { return 0.5 * kLambda * U_exact(x, y); }
inline double F_rhs(double x, double y) { return 0.5 * kLambda * kLambda * U_exact(x, y); }
inline arma::uword idx(u32 i, u32 j, u32 m) { return j * (m + 2) + i; }
}

using arma::sp_mat;
using arma::vec;

static void fill_tensor_identity_2D(u32 m, u32 n, TensorField<2>& K) {
  const arma::uword N = (m + 2) * (n + 2);
  K.xx.set_size(N); K.yy.set_size(N);
  K.xy.set_size(N); K.yx.set_size(N);
  K.xx.ones(); K.yy.ones();
  K.xy.zeros(); K.yx.zeros();
}

int main() {
  const std::vector<u32> Mvals = {10, 20, 40, 80};
  constexpr u16 kOrder = 2;

  const double alpha = -std::exp(kLambda);
  const double beta  = (std::exp(kLambda) - 1.0) / kLambda;

  std::cout << std::setprecision(10) << std::scientific;
  std::cout << "Grid Size\tL2 Rel Error\tRate\n";

  double prev_err = 0.0;

  for (size_t r = 0; r < Mvals.size(); ++r) {
    const u32 m = Mvals[r], n = Mvals[r];
    const double dx = Lx / double(m);
    const double dy = Ly / double(n);

    arma::vec Xcc(m + 2), Ycc(n + 2);
    Xcc(0) = -0.5 * dx;            Xcc(m + 1) = 1.0 + 0.5 * dx;
    for (u32 i = 1; i <= m; ++i) Xcc(i) = (i - 0.5) * dx;
    Ycc(0) = -0.5 * dy;            Ycc(n + 1) = 1.0 + 0.5 * dy;
    for (u32 j = 1; j <= n; ++j) Ycc(j) = (j - 0.5) * dy;

    const arma::uword N = (m + 2) * (n + 2);

    TensorField<2> K;
    fill_tensor_identity_2D(m, n, K);

    FluxND<2> flux(K);
    sp_mat A = flux.diffusion(kOrder, m, n, dx, dy);

    sp_mat B = static_cast<sp_mat>(RobinBC(kOrder, m, dx, n, dy, alpha, beta));
    sp_mat System = A + B;

    vec rhs(N, arma::fill::zeros);
    for (u32 j = 1; j <= n; ++j)
      for (u32 i = 1; i <= m; ++i)
        rhs(idx(i, j, m)) = F_rhs(Xcc(i), Ycc(j));

    vec gvec(N, arma::fill::zeros);
    for (u32 j = 1; j <= n; ++j) {
      const double ue0 = U_exact(0.0, Ycc(j)), ux0 = Ux_exact(0.0, Ycc(j));
      const double ue1 = U_exact(1.0, Ycc(j)), ux1 = Ux_exact(1.0, Ycc(j));
      gvec(idx(0, j, m))     += alpha * ue0 + beta * (-ux0);
      gvec(idx(m + 1, j, m)) += alpha * ue1 + beta * (+ux1);
    }
    for (u32 i = 1; i <= m; ++i) {
      const double ue0 = U_exact(Xcc(i), 0.0), uy0 = Uy_exact(Xcc(i), 0.0);
      const double ue1 = U_exact(Xcc(i), 1.0), uy1 = Uy_exact(Xcc(i), 1.0);
      gvec(idx(i, 0, m))     += alpha * ue0 + beta * (-uy0);
      gvec(idx(i, n + 1, m)) += alpha * ue1 + beta * (+uy1);
    }

    rhs += gvec;

    vec U;
    bool ok = arma::spsolve(U, System, rhs, "superlu");
    if (!ok) U = arma::spsolve(System, rhs);

    double err2 = 0.0, ref2 = 0.0;
    for (u32 j = 1; j <= n; ++j)
      for (u32 i = 1; i <= m; ++i) {
        const double ue = U_exact(Xcc(i), Ycc(j));
        const double du = U(idx(i, j, m)) - ue;
        err2 += du * du;
        ref2 += ue * ue;
      }
    const double relL2 = std::sqrt(err2 * dx * dy) / std::sqrt(ref2 * dx * dy);

    std::cout << std::setw(2) << m << "x" << n << "\t" << relL2;
    if (r > 0) std::cout << "\t" << std::log2(prev_err / relL2);
    std::cout << "\n";

    prev_err = relL2;
  }
  return 0;
}
