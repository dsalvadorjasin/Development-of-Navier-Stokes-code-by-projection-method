// C++ counterpart of variables_module.f90: the global simulation state.
//
// Compile-time parameters are constexpr; everything else is a global that
// keeps the exact Fortran name and index range.  Scalars are exposed as
// references and arrays as FArray views so that the same routine code works
// in both builds:
//   * pure C++ build   : storage owned by variables.cpp
//   * hybrid build     : views over the gfortran module storage
//                        (-DNS_HYBRID, see cpp.mk / variables.cpp)
//
// Fortran declared xx, yy as bare `real`; the Makefile compiles with
// -fdefault-real-8 so they are double precision too.
#pragma once

#include "fortran_array.hpp"

namespace ns {

constexpr int nx = 11, ny = 11, nz = 1;
constexpr double nu = 1e-02, dt2 = 1e-04;
constexpr int ndim = nx * ny, mdim = 3;

// integers
extern int &itmax, &isto, &istep, &nstep, &Q;
extern FArray1<int> jcoef;  // (1:mdim)
extern FArray1<int> Ldiag;  // (1:5)

// reals
extern double &max_x, &max_y, &min_u, &min_v;
extern double &zeta, &time, &dt, &dx, &dy, &ddx, &ddt, &ddy, &nu1;

// fields with ghost cells (0:nx+1,0:ny+1)
extern FArray2<double> p, u, v, u1, v1, u_star, v_star, abs_u, abs_v;
// cell-centred fields (1:nx,1:ny)
extern FArray2<double> div, uc, vc, rot, w, pre;
// sparse matrix in ELL/diagonal storage (1:ndim,1:mdim) and its LU factor (1:ndim,1:5)
extern FArray2<double> coef, L;
// solver work vectors (1:ndim)
extern FArray1<double> div1, p_s, r_s, r2_s, q_s, s_s, x1;
// cell-centre coordinates (1:nx), (1:ny)
extern FArray1<double> xx, yy;

}  // namespace ns
