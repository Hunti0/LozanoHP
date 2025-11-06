
#include <armadillo>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "mole.h"
#include "Flux3D.h"   

using arma::sp_mat;
using arma::vec;
using arma::uword;

using u16 = std::uint16_t;

// ---------------- Grid & indexing ----------------
static inline uword idx(unsigned i, unsigned j, unsigned m) {
  return i + (m + 2) * j;   
}
struct Grid {
  unsigned m, n;
  double dx, dy;
  std::vector<double> x, y; 
};
static Grid make_grid(unsigned m, unsigned n) {
  Grid g{m, n, 1.0/m, 1.0/n, {}, {}};
  g.x.resize(m + 2); g.y.resize(n + 2);
  g.x[0] = 0.0;      g.x[m+1] = 1.0;
  g.y[0] = 0.0;      g.y[n+1] = 1.0;
  for (unsigned i = 1; i <= m; ++i) g.x[i] = (i - 0.5) * g.dx;
  for (unsigned j = 1; j <= n; ++j) g.y[j] = (j - 0.5) * g.dy;
  return g;
}

// ---------------- Problem data ----------------
static inline double u_exact(double x, double y) { return std::exp(x * y); }

struct Kconst { double a,b,c,d; };
static inline double f_rhs(double x, double y, const Kconst& K) {
  const double exy = std::exp(x*y);
  return -(K.d*x*x + K.a*y*y + (K.b+K.c)*(1.0 + x*y)) * exy;
}

// ---------------- Assembly: A = -D F G ----------------
static sp_mat assemble_A(u16 k, unsigned m, unsigned n, double dx, double dy,
                         const TensorField<2>& KT)
{
  FluxND<2> flux(KT);
  sp_mat F  = flux.flux_matrix(k, m, n);
  sp_mat Gm = Gradient(k, m, n, dx, dy);
  sp_mat Dm = Divergence(k, m, n, dx, dy);
  return -(Dm * F * Gm);
}

// ---------------- Boundary tooling ----------------

// Build the boundary index list (all nodes with i∈{0,m+1} or j∈{0,n+1})
static std::vector<uword> boundary_rows(unsigned m, unsigned n) {
  std::vector<uword> B; B.reserve(2*(m+n+2));
  for (unsigned j = 0; j <= n + 1; ++j) {
    B.push_back(idx(0,     j, m));
    B.push_back(idx(m + 1, j, m));
  }
  for (unsigned i = 0; i <= m + 1; ++i) {
    B.push_back(idx(i, 0,     m));
    B.push_back(idx(i, n + 1, m));
  }
  return B;
}

static sp_mat build_normal_picker(unsigned m, unsigned n) {
  const arma::uword Nc = (m + 2) * (n + 2);
  const arma::uword Nu = (m + 1) * n;      
  const arma::uword Nv = m * (n + 1);     

  auto u_face = [m](unsigned i, unsigned j)->arma::uword { return j * (m + 1) + i; };
  auto v_face = [m, n](unsigned i, unsigned j)->arma::uword { return i * (n + 1) + j; };

  sp_mat N(Nc, Nu + Nv);

  
  for (unsigned j = 1; j <= n; ++j) {
    arma::uword row = idx(0, j, m);
    arma::uword uf  = u_face(0, j - 1);
    N(row, uf) = -1.0;    // n = (-1,0)
  }
  
  for (unsigned j = 1; j <= n; ++j) {
    arma::uword row = idx(m + 1, j, m);
    arma::uword uf  = u_face(m, j - 1);
    N(row, uf) = +1.0;    // n = (+1,0)
  }

  for (unsigned i = 1; i <= m; ++i) {
    arma::uword row = idx(i, 0, m);
    arma::uword vf  = v_face(i - 1, 0);
    N(row, Nu + vf) = -1.0; // n = (0,-1)
  }
  
  for (unsigned i = 1; i <= m; ++i) {
    arma::uword row = idx(i, n + 1, m);
    arma::uword vf  = v_face(i - 1, n);
    N(row, Nu + vf) = +1.0; // n = (0,+1)
  }
  return N;
}

static sp_mat build_boundary_row_scaler(unsigned m, unsigned n,
                                        const TensorField<2>& KT) {
  const arma::uword Nc = (m + 2) * (n + 2);
  vec s(Nc, arma::fill::zeros);

 
  for (unsigned j = 1; j <= n; ++j) {
    s(idx(0,     j, m)) = KT.xx(idx(0,     j, m));
    s(idx(m + 1, j, m)) = KT.xx(idx(m + 1, j, m));
  }
  
  for (unsigned i = 1; i <= m; ++i) {
    s(idx(i, 0,     m)) = KT.yy(idx(i, 0,     m));
    s(idx(i, n + 1, m)) = KT.yy(idx(i, n + 1, m));
  }

  
  sp_mat S(Nc, Nc);
  for (arma::uword p = 0; p < Nc; ++p) if (s(p) != 0.0) S(p,p) = s(p);
  return S;
}

static sp_mat build_robin_enhanced(u16 k, unsigned m, unsigned n, double dx, double dy,
                                   const TensorField<2>& KT, double a, double b)
{
  
  RobinBC RBC_id(k, m, dx, n, dy, 1.0, 0.0); 
  RobinBC RBC_dn(k, m, dx, n, dy, 0.0, 1.0); 


  sp_mat S = build_boundary_row_scaler(m, n, KT);
  sp_mat R_norm = a * sp_mat(RBC_id) + b * (S * sp_mat(RBC_dn));


  sp_mat G = Gradient(k, m, n, dx, dy);

  // Full flux with full K
  FluxND<2> flux_full(KT);
  sp_mat F_full = flux_full.flux_matrix(k, m, n);

  TensorField<2> Kdiag{KT.xx, vec(KT.xy.n_rows, arma::fill::zeros),
                       vec(KT.yx.n_rows, arma::fill::zeros), KT.yy};
  FluxND<2> flux_diag(Kdiag);
  sp_mat F_diag = flux_diag.flux_matrix(k, m, n);

  sp_mat F_off = F_full - F_diag;   

  sp_mat Npick = build_normal_picker(m, n);
  sp_mat R_off = b * (Npick * F_off * G);

  return R_norm + R_off;
}

static inline void apply_bc_rows_mat(const sp_mat& A0, const sp_mat& R,
                                     const vec& rhs_interior, const vec& uex,
                                     unsigned m, unsigned n,
                                     sp_mat& A_out, vec& b_out)
{
  const uword Nc  = A0.n_rows;

  sp_mat Aint = A0;
  const auto B = boundary_rows(m, n);
  for (uword p : B) Aint.row(p).zeros();

  A_out = Aint + R;

  b_out = rhs_interior;
  vec g = R * uex;
  for (uword p : B) b_out(p) = g(p);

  if (A_out.n_rows != Nc || A_out.n_cols != Nc || b_out.n_rows != Nc)
    throw std::runtime_error("apply_bc_rows_mat: size mismatch");
}


static inline std::pair<double,double> l2_and_linf(const vec& u,
                                                   const vec& ue,
                                                   unsigned m, unsigned n,
                                                   double, double)
{
  const double Nc = double((m + 2) * (n + 2));
  double L2   = arma::norm(u - ue, 2) / std::sqrt(Nc);
  double Linf = arma::norm(u - ue, "inf");
  return {L2, Linf};
}

static inline double observed_order(double E1, double E2, double h1, double h2) {
  return std::log(E1 / E2) / std::log(h1 / h2);
}

// =========================== main ===========================
int main() try {
  const u16 kOrder = 2;
  const Kconst K{2.0, 1.0, 1.0, 2.0};   // [[2,1],[1,2]]
  std::vector<unsigned> Ms = {10, 17, 20, 33, 65};

  std::vector<double> L2_dir, Linf_dir, L2_robin, Linf_robin;

  for (unsigned m : Ms) {
    unsigned n = m;
    Grid g = make_grid(m, n);
    const uword Nc = (m + 2) * (n + 2);

    // tensor fields 
    vec Kxx(Nc, arma::fill::value(K.a));
    vec Kxy(Nc, arma::fill::value(K.b));
    vec Kyx(Nc, arma::fill::value(K.c));
    vec Kyy(Nc, arma::fill::value(K.d));
    TensorField<2> KT{Kxx, Kxy, Kyx, Kyy};

    // assemble interior operator and interior RHS f
    sp_mat A0 = assemble_A(kOrder, m, n, g.dx, g.dy, KT);
    vec rhs_interior(Nc, arma::fill::zeros);
    for (unsigned j = 1; j <= n; ++j)
      for (unsigned i = 1; i <= m; ++i)
        rhs_interior(idx(i,j,m)) = f_rhs(g.x[i], g.y[j], K);


    vec uex(Nc, arma::fill::zeros);
    for (unsigned j = 0; j <= n + 1; ++j)
      for (unsigned i = 0; i <= m + 1; ++i)
        uex(idx(i,j,m)) = u_exact(g.x[i], g.y[j]);

    // --- Robin ---
    {
      sp_mat R_eff = build_robin_enhanced(kOrder, m, n, g.dx, g.dy, KT, /*a=*/1.0, /*b=*/1.0);
      sp_mat A; vec b;
      apply_bc_rows_mat(A0, R_eff, rhs_interior, uex, m, n, A, b);
      vec u; bool ok = arma::spsolve(u, A, b, "superlu");
      if (!ok) { std::cerr << "Robin solve failed at n=" << m << "\n"; return 1; }
      auto [L2, Linf] = l2_and_linf(u, uex, m, n, g.dx, g.dy);
      L2_robin.push_back(L2); Linf_robin.push_back(Linf);
    }

    // --- Dirichlet ---
    {
      RobinBC RBC_dir(kOrder, m, g.dx, n, g.dy, 1.0, 0.0);
      sp_mat R_dir = static_cast<sp_mat>(RBC_dir);
      sp_mat A; vec b;
      apply_bc_rows_mat(A0, R_dir, rhs_interior, uex, m, n, A, b);
      vec u; bool ok = arma::spsolve(u, A, b, "superlu");
      if (!ok) { std::cerr << "Dirichlet solve failed at n=" << m << "\n"; return 1; }
      auto [L2, Linf] = l2_and_linf(u, uex, m, n, g.dx, g.dy);
      L2_dir.push_back(L2); Linf_dir.push_back(Linf);
    }
  }


  std::vector<unsigned> Ns = Ms;
  std::vector<double> hs; hs.reserve(Ns.size());
  for (auto n : Ns) hs.push_back(1.0 / double(n));

  auto fmtE = [](double x){
    std::ostringstream ss; ss<<std::scientific<<std::uppercase<<std::setprecision(2)<<x;
    return ss.str();
  };

  auto print_table = [&](const char* title,
                         const std::vector<double>& L2v,
                         const std::vector<double>& Linfv) {
    using std::cout;
    cout << "\n" << title << "\n";
    cout << std::left << std::setw(6) << "n"
         << std::right << std::setw(14) << "L2"
         << std::setw(10) << "Order"
         << std::setw(14) << "Max"
         << std::setw(10) << "Order" << "\n";
    for (size_t i=0;i<Ns.size();++i) {
      double p_L2   = (i==0 ? NAN : observed_order(L2v[i-1],  L2v[i],  hs[i-1], hs[i]));
      double p_Linf = (i==0 ? NAN : observed_order(Linfv[i-1],Linfv[i],hs[i-1], hs[i]));
      std::ostringstream s1,s2;
      if (i>0) { s1<<std::fixed<<std::setprecision(4)<<p_L2; s2<<std::fixed<<std::setprecision(4)<<p_Linf; }
      cout << std::left << std::setw(6) << Ns[i]
           << std::right << std::setw(14) << fmtE(L2v[i])
           << std::setw(10) << (i==0 ? "" : s1.str())
           << std::setw(14) << fmtE(Linfv[i])
           << std::setw(10) << (i==0 ? "" : s2.str())
           << "\n";
    }
  };

  print_table("Results for Problem 4 — Robin", L2_robin, Linf_robin);
  print_table("Results for Problem 4 — Dirichlet", L2_dir, Linf_dir);

  return 0;
}
catch (const std::exception& e) {
  std::cerr << "ERROR: " << e.what() << "\n";
  return 1;
}

