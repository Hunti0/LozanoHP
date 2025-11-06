#ifndef FLUX3D_H
#define FLUX3D_H

#include <stdexcept>
#include <string>
#include <vector>
#include <armadillo>

#include "utils.h"
#include "mole.h"

using arma::sp_mat;
using arma::vec;
using arma::uword;

using u16 = std::uint16_t;
using u32 = std::uint32_t;

//======================== Tensor fields ========================
template<int Dim> struct TensorField;

template<> struct TensorField<1> {
  vec k;
  static TensorField<1> Scalar(const vec& k_)    { return {k_}; }
  static TensorField<1> Isotropic(const vec& k_) { return {k_}; }
};

template<> struct TensorField<2> {
  vec xx, xy, yx, yy;
  static TensorField<2> Isotropic(const vec& k) {
    vec z(k.n_elem, arma::fill::zeros);
    return {k, z, z, k};
  }
  static TensorField<2> Diagonal(const vec& kxx, const vec& kyy) {
    uword N = kxx.n_elem; vec z(N, arma::fill::zeros);
    return {kxx, z, z, kyy};
  }
};

//======================== Small helpers ========================
inline vec elem_recip(const vec& x) {
  vec r = x;
  r.transform([](double a) {
    const double eps = 1e-14;
    return 1.0 / (std::abs(a) < eps ? (a >= 0 ? eps : -eps) : a);
  });
  return r;
}
inline sp_mat Diag(const vec& v) { sp_mat D(v.n_rows, v.n_rows); D.diag() = v; return D; }
inline void ensure_len(const vec& v, uword expect, const char* name) {
  if (v.n_elem != expect) {
    throw std::invalid_argument(std::string("Tensor component '")+name+
      "' has length " + std::to_string(v.n_elem) + ", expected " + std::to_string(expect) + ".");
  }
}

//======================== Face blocks (2D) ========================
struct FaceBlocks2D {
  sp_mat CU, CV;
  sp_mat UC, VC;
  sp_mat VU, UV;
};

inline sp_mat CF_1D(u32 m) {
  return Interpol(m, 0.5);
}

inline FaceBlocks2D blocks_2D(u32 m, u32 n) {
  Interpol CF(m, n, 0.5, 0.5);
  uword U = (m + 1) * n;
  uword V = m * (n + 1);
  sp_mat CU = CF.rows(0, U - 1);
  sp_mat CV = CF.rows(U, U + V - 1);

  Interpol FC(true, m, n, 0.5, 0.5);
  sp_mat UC = FC.cols(0, U - 1);
  sp_mat VC = FC.cols(U, U + V - 1);

  sp_mat VU = CU * VC;
  sp_mat UV = CV * UC;
  return {CU, CV, UC, VC, VU, UV};
}

//======================== Tangential face↔face maps (built here) ========================
inline sp_mat Ax_interp(uword m, u16 k) {
  const uword rows = m + 1, cols = m;
  sp_mat A(rows, cols);
  if (m == 0) return A;

  if (k < 4 || m < 4) {
    A(0,0) = 1.0;
    for (uword i=1;i<=m-1;++i){ A(i,i-1)=0.5; A(i,i)=0.5; }
    A(m,m-1) = 1.0;
    return A;
  }

  A(0,0) =  35.0/16.0;  A(0,1) = -35.0/16.0;
  A(0,2) =  21.0/16.0;  A(0,3) =  -5.0/16.0;

  A(1,0) =   5.0/16.0;  A(1,1) =  15.0/16.0;
  A(1,2) =  -5.0/16.0;  A(1,3) =   1.0/16.0;

  for (uword i=2; i<=m-2; ++i) {
    A(i,i-2) = -1.0/16.0;  A(i,i-1) =  9.0/16.0;
    A(i,i)   =  9.0/16.0;  A(i,i+1) = -1.0/16.0;
  }

  A(m-1, m-4) =  -5.0/16.0;  A(m-1, m-3) =  21.0/16.0;
  A(m-1, m-2) = -35.0/16.0;  A(m-1, m-1) =  35.0/16.0;

  A(m,   m-4) =  -5.0/16.0;  A(m,   m-3) =  21.0/16.0;
  A(m,   m-2) = -35.0/16.0;  A(m,   m-1) =  35.0/16.0;

  return A;
}

inline sp_mat Ay_interp(uword n, u16 k) {
  const uword rows = n + 1, cols = n;
  sp_mat A(rows, cols);
  if (n == 0) return A;

  if (k < 4 || n < 4) {
    A(0,0) = 1.0;
    for (uword j=1;j<=n-1;++j){ A(j,j-1)=0.5; A(j,j)=0.5; }
    A(n,n-1) = 1.0;
    return A;
  }

  A(0,0) =  35.0/16.0;  A(0,1) = -35.0/16.0;
  A(0,2) =  21.0/16.0;  A(0,3) =  -5.0/16.0;

  A(1,0) =   5.0/16.0;  A(1,1) =  15.0/16.0;
  A(1,2) =  -5.0/16.0;  A(1,3) =   1.0/16.0;

  for (uword j=2; j<=n-2; ++j) {
    A(j,j-2) = -1.0/16.0;  A(j,j-1) =  9.0/16.0;
    A(j,j)   =  9.0/16.0;  A(j,j+1) = -1.0/16.0;
  }

  A(n-1, n-4) =  -5.0/16.0;  A(n-1, n-3) =  21.0/16.0;
  A(n-1, n-2) = -35.0/16.0;  A(n-1, n-1) =  35.0/16.0;

  A(n,   n-4) =  -5.0/16.0;  A(n,   n-3) =  21.0/16.0;
  A(n,   n-2) = -35.0/16.0;  A(n,   n-1) =  35.0/16.0;

  return A;
}

inline sp_mat Ix_VtoU(u16 k, uword m, uword n) {
  sp_mat Axp = Ax_interp(m, k);
  sp_mat Aym = Ay_interp(n, k).t();
  return arma::kron(Aym, Axp);
}
inline sp_mat Iy_UtoV(u16 k, uword m, uword n) {
  sp_mat Axm = Ax_interp(m, k).t();
  sp_mat Ayp = Ay_interp(n, k);
  return arma::kron(Ayp, Axm);
}

//======================== FluxND ========================
template<int Dim> class FluxND;

// ---------- 1D ----------
template<> class FluxND<1> {
public:
  explicit FluxND(const TensorField<1>& K) : K_(K) {}
  sp_mat flux_matrix(u16 /*k*/, u32 m) const {
    const uword Nc = (m + 2);
    ensure_len(K_.k, Nc, "k");
    sp_mat CF = CF_1D(m);
    vec inv_k  = elem_recip(K_.k);
    vec inv_kf = CF * inv_k;
    vec kf     = elem_recip(inv_kf);
    return Diag(kf);
  }
  sp_mat diffusion(u16 k, u32 m, Real dx) const {
    Gradient   G(k, m, dx);
    Divergence Dv(k, m, dx);
    return (sp_mat)Dv * flux_matrix(k, m) * (sp_mat)G;
  }
private:
  TensorField<1> K_;
};

// ---------- 2D ----------
template<> class FluxND<2> {
public:
  explicit FluxND(const TensorField<2>& K) : K_(K) {}
  sp_mat flux_matrix(u16 k, u32 m, u32 n) const {
    const uword Nc = (m + 2) * (n + 2);
    ensure_len(K_.xx, Nc, "xx"); ensure_len(K_.xy, Nc, "xy");
    ensure_len(K_.yx, Nc, "yx"); ensure_len(K_.yy, Nc, "yy");

    FaceBlocks2D B = blocks_2D(m, n);
    vec xxU = elem_recip( B.CU * elem_recip(K_.xx) );
    vec yyV = elem_recip( B.CV * elem_recip(K_.yy) );
    vec xyU = B.CU * K_.xy;
    vec yxV = B.CV * K_.yx;

    sp_mat Ix = Ix_VtoU(k, m, n);
    sp_mat Iy = Iy_UtoV(k, m, n);

    sp_mat Fxx = Diag(xxU);
    sp_mat Fyy = Diag(yyV);
    sp_mat Fxy = Diag(xyU) * Ix;
    sp_mat Fyx = Diag(yxV) * Iy;

    return arma::join_rows( arma::join_cols(Fxx, Fyx),
                            arma::join_cols(Fxy, Fyy) );
  }
  sp_mat diffusion(u16 k, u32 m, u32 n, Real dx, Real dy) const {
    Gradient   G(k, m, n, dx, dy);
    Divergence Dv(k, m, n, dx, dy);
    return (sp_mat)Dv * flux_matrix(k, m, n) * (sp_mat)G;
  }
private:
  TensorField<2> K_;
};

#endif
