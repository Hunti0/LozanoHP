#include <armadillo>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <algorithm>

#include "Flux3D.h"
#include "robinbc.h"

using arma::sp_mat;
using arma::vec;
using std::cout;
using u16 = unsigned short;
using u32 = unsigned int;
using Real = double;

inline arma::uword idx(u32 i, u32 j, u32 m_plus_2) { return i + (arma::uword)m_plus_2 * j; }

static std::vector<Real> cwb_coords(u32 m, Real L = 1.0) {
  std::vector<Real> x(m + 2);
  Real dx = L / static_cast<Real>(m);
  x[0] = 0.0;
  for (u32 i = 1; i <= m; ++i) x[i] = (i - 0.5) * dx;
  x[m + 1] = L;
  return x;
}

inline Real u_exact(Real x, Real y) {
  return std::pow(x, 3) * y + std::pow(y, 4) + std::sin(x) * std::cos(y);
}
inline Real divKgrad_u(Real x, Real y) {
  return 60.0 * x * y + 12.0 * y * y - 11.0 * std::sin(x) * std::cos(y);
}
inline Real f_paper(Real x, Real y) { return -divKgrad_u(x, y); }

static TensorField<2> makeKdiag(u32 m, u32 n, Real Kxx = 10.0, Real Kyy = 1.0) {
  arma::uword Nc = (m + 2) * (n + 2);
  vec kxx(Nc), kyy(Nc);
  kxx.fill(Kxx);
  kyy.fill(Kyy);
  return TensorField<2>::Diagonal(kxx, kyy);
}

static vec build_vec_by_center(u32 m, u32 n,
                               const std::vector<Real>& X,
                               const std::vector<Real>& Y,
                               Real (*fun)(Real, Real)) {
  vec v((m + 2) * (n + 2), arma::fill::zeros);
  for (u32 j = 0; j <= n + 1; ++j)
    for (u32 i = 0; i <= m + 1; ++i)
      v(idx(i, j, m + 2)) = fun(X[i], Y[j]);
  return v;
}

static vec build_vec_interior(u32 m, u32 n,
                              const std::vector<Real>& X,
                              const std::vector<Real>& Y,
                              Real (*fun)(Real, Real)) {
  vec v((m + 2) * (n + 2), arma::fill::zeros);
  for (u32 j = 1; j <= n; ++j)
    for (u32 i = 1; i <= m; ++i)
      v(idx(i, j, m + 2)) = fun(X[i], Y[j]);
  return v;
}

static void apply_bc_rows(const sp_mat& A0, const sp_mat& RBC,
                          const vec& interior_rhs, const vec& uex,
                          sp_mat& A, vec& b) {
  const arma::uword N = A0.n_rows;
  arma::Col<unsigned char> has_nz(N, arma::fill::zeros);
  for (sp_mat::const_iterator it = RBC.begin(); it != RBC.end(); ++it) has_nz(it.row()) = 1;

  sp_mat I = arma::speye(N, N), P(N, N);
  for (arma::uword r = 0; r < N; ++r) if (has_nz(r)) P(r, r) = 1.0;

  A = (I - P) * A0 + RBC;
  b = (I - P) * interior_rhs + RBC * uex;
}

static std::pair<double, double>
l2_and_linf(const vec& u, const vec& uex, u32 m, u32 n, double dx, double dy) {
  double sum = 0.0, linf = 0.0;
  for (u32 j = 1; j <= n; ++j)
    for (u32 i = 1; i <= m; ++i) {
      double e = u(idx(i, j, m + 2)) - uex(idx(i, j, m + 2));
      sum += e * e;
      linf = std::max(linf, std::abs(e));
    }
  return {std::sqrt(sum * dx * dy), linf};
}

static double observed_order(double e_coarse, double e_fine, double h_coarse, double h_fine) {
  return std::log(e_coarse / e_fine) / std::log(h_coarse / h_fine);
}

int main() {
  arma::arma_rng::set_seed(1);
  const u16 kOrder = 2;
  const std::vector<u32> Ns = {10, 17, 20, 33, 65};
  const double L = 1.0;

  std::vector<double> hs, L2_robin, Linf_robin, L2_dir, Linf_dir;

  for (size_t t = 0; t < Ns.size(); ++t) {
    u32 m = Ns[t], n = Ns[t];
    Real dx = L / static_cast<Real>(m), dy = L / static_cast<Real>(n);
    hs.push_back(std::max(dx, dy));

    auto X = cwb_coords(m, L);
    auto Y = cwb_coords(n, L);
    auto K = makeKdiag(m, n, 10.0, 1.0);

    FluxND<2> flux(K);
    sp_mat A0 = flux.diffusion(kOrder, m, n, dx, dy);

    vec uex = build_vec_by_center(m, n, X, Y, u_exact);
    vec f_pap = build_vec_interior(m, n, X, Y, f_paper);
    vec rhs_interior = -f_pap;

    {
      RobinBC RBC_robin(kOrder, m, dx, n, dy, 1.0, 1.0);
      sp_mat A, bRBC = static_cast<sp_mat>(RBC_robin);
      vec b;
      apply_bc_rows(A0, bRBC, rhs_interior, uex, A, b);
      vec u; bool ok = arma::spsolve(u, A, b, "superlu");
      if (!ok) { std::cerr << "Robin solve failed at n=" << m << "\n"; return 1; }
      auto [L2, Linf] = l2_and_linf(u, uex, m, n, dx, dy);
      L2_robin.push_back(L2); Linf_robin.push_back(Linf);
    }
    {
      RobinBC RBC_dir(kOrder, m, dx, n, dy, 1.0, 0.0);
      sp_mat A, bRBC = static_cast<sp_mat>(RBC_dir);
      vec b;
      apply_bc_rows(A0, bRBC, rhs_interior, uex, A, b);
      vec u; bool ok = arma::spsolve(u, A, b, "superlu");
      if (!ok) { std::cerr << "Dirichlet solve failed at n=" << m << "\n"; return 1; }
      auto [L2, Linf] = l2_and_linf(u, uex, m, n, dx, dy);
      L2_dir.push_back(L2); Linf_dir.push_back(Linf);
    }
  }

  auto fmtE = [](double x){ std::ostringstream ss; ss<<std::scientific<<std::uppercase<<std::setprecision(2)<<x; return ss.str(); };
  auto print_table = [&](const char* title, const std::vector<double>& L2v, const std::vector<double>& Linfv) {
    cout << "\n" << title << "\n";
    cout << std::left << std::setw(6) << "n"
         << std::right << std::setw(14) << "L2"
         << std::setw(10) << "Order"
         << std::setw(14) << "Max"
         << std::setw(10) << "Order" << "\n";
    for (size_t i = 0; i < Ns.size(); ++i) {
      double pL2   = (i == 0 ? NAN : observed_order(L2v[i-1],  L2v[i],  hs[i-1], hs[i]));
      double pLinf = (i == 0 ? NAN : observed_order(Linfv[i-1],Linfv[i],hs[i-1], hs[i]));
      std::ostringstream s1, s2;
      if (i > 0) { s1 << std::fixed << std::setprecision(4) << pL2; s2 << std::fixed << std::setprecision(4) << pLinf; }
      cout << std::left << std::setw(6) << Ns[i]
           << std::right << std::setw(14) << fmtE(L2v[i])
           << std::setw(10) << (i == 0 ? "" : s1.str())
           << std::setw(14) << fmtE(Linfv[i])
           << std::setw(10) << (i == 0 ? "" : s2.str())
           << "\n";
    }
  };

  print_table("Results for Problem 3 — Robin (a=b=1)", L2_robin, Linf_robin);
  print_table("Results for Problem 3 — Dirichlet (a=1,b=0)", L2_dir, Linf_dir);
  return 0;
}
