// Storage for the globals declared in variables.hpp.
//
// Pure C++ build : static storage owned here (zero-initialised like Fortran
//                  BSS module variables).
// Hybrid build   : -DNS_HYBRID makes every global an alias of the gfortran
//                  module variable (mangled `__variables_MOD_<name>`), so
//                  C++ routines linked into the Fortran program read and
//                  write exactly the same memory as the remaining Fortran
//                  routines.
#include "variables.hpp"

namespace ns {

#ifdef NS_HYBRID

extern "C" {
extern int __variables_MOD_itmax, __variables_MOD_isto, __variables_MOD_istep,
    __variables_MOD_nstep, __variables_MOD_q;
extern int __variables_MOD_jcoef[mdim], __variables_MOD_ldiag[5];
extern double __variables_MOD_max_x, __variables_MOD_max_y, __variables_MOD_min_u,
    __variables_MOD_min_v, __variables_MOD_zeta, __variables_MOD_time, __variables_MOD_dt,
    __variables_MOD_dx, __variables_MOD_dy, __variables_MOD_ddx, __variables_MOD_ddt,
    __variables_MOD_ddy, __variables_MOD_nu1;
extern double __variables_MOD_p[], __variables_MOD_u[], __variables_MOD_v[],
    __variables_MOD_u1[], __variables_MOD_v1[], __variables_MOD_u_star[],
    __variables_MOD_v_star[], __variables_MOD_abs_u[], __variables_MOD_abs_v[];
extern double __variables_MOD_div[], __variables_MOD_uc[], __variables_MOD_vc[],
    __variables_MOD_rot[], __variables_MOD_w[], __variables_MOD_pre[];
extern double __variables_MOD_coef[], __variables_MOD_l[];
extern double __variables_MOD_div1[], __variables_MOD_p_s[], __variables_MOD_r_s[],
    __variables_MOD_r2_s[], __variables_MOD_q_s[], __variables_MOD_s_s[], __variables_MOD_x1[];
extern double __variables_MOD_xx[], __variables_MOD_yy[];
}
#define NS_SYM(name) __variables_MOD_##name
#define NS_SYM_Q __variables_MOD_q
#define NS_SYM_L __variables_MOD_l

#else

namespace {
constexpr int ng = (nx + 2) * (ny + 2);
int itmax_store, isto_store, istep_store, nstep_store, Q_store;
int jcoef_store[mdim], ldiag_store[5];
double max_x_store, max_y_store, min_u_store, min_v_store, zeta_store, time_store, dt_store,
    dx_store, dy_store, ddx_store, ddt_store, ddy_store, nu1_store;
double p_store[ng], u_store[ng], v_store[ng], u1_store[ng], v1_store[ng], u_star_store[ng],
    v_star_store[ng], abs_u_store[ng], abs_v_store[ng];
double div_store[nx * ny], uc_store[nx * ny], vc_store[nx * ny], rot_store[nx * ny],
    w_store[nx * ny], pre_store[nx * ny];
double coef_store[ndim * mdim], L_store[ndim * 5];
double div1_store[ndim], p_s_store[ndim], r_s_store[ndim], r2_s_store[ndim], q_s_store[ndim],
    s_s_store[ndim], x1_store[ndim];
double xx_store[nx], yy_store[ny];
}  // namespace
#define NS_SYM(name) name##_store
#define NS_SYM_Q Q_store
#define NS_SYM_L L_store

#endif

int &itmax = NS_SYM(itmax), &isto = NS_SYM(isto), &istep = NS_SYM(istep), &nstep = NS_SYM(nstep),
    &Q = NS_SYM_Q;
FArray1<int> jcoef(NS_SYM(jcoef), 1, mdim);
FArray1<int> Ldiag(NS_SYM(ldiag), 1, 5);

double &max_x = NS_SYM(max_x), &max_y = NS_SYM(max_y), &min_u = NS_SYM(min_u),
       &min_v = NS_SYM(min_v);
double &zeta = NS_SYM(zeta), &time = NS_SYM(time), &dt = NS_SYM(dt), &dx = NS_SYM(dx),
       &dy = NS_SYM(dy), &ddx = NS_SYM(ddx), &ddt = NS_SYM(ddt), &ddy = NS_SYM(ddy),
       &nu1 = NS_SYM(nu1);

FArray2<double> p(NS_SYM(p), 0, nx + 1, 0, ny + 1);
FArray2<double> u(NS_SYM(u), 0, nx + 1, 0, ny + 1);
FArray2<double> v(NS_SYM(v), 0, nx + 1, 0, ny + 1);
FArray2<double> u1(NS_SYM(u1), 0, nx + 1, 0, ny + 1);
FArray2<double> v1(NS_SYM(v1), 0, nx + 1, 0, ny + 1);
FArray2<double> u_star(NS_SYM(u_star), 0, nx + 1, 0, ny + 1);
FArray2<double> v_star(NS_SYM(v_star), 0, nx + 1, 0, ny + 1);
FArray2<double> abs_u(NS_SYM(abs_u), 0, nx + 1, 0, ny + 1);
FArray2<double> abs_v(NS_SYM(abs_v), 0, nx + 1, 0, ny + 1);

FArray2<double> div(NS_SYM(div), 1, nx, 1, ny);
FArray2<double> uc(NS_SYM(uc), 1, nx, 1, ny);
FArray2<double> vc(NS_SYM(vc), 1, nx, 1, ny);
FArray2<double> rot(NS_SYM(rot), 1, nx, 1, ny);
FArray2<double> w(NS_SYM(w), 1, nx, 1, ny);
FArray2<double> pre(NS_SYM(pre), 1, nx, 1, ny);

FArray2<double> coef(NS_SYM(coef), 1, ndim, 1, mdim);
FArray2<double> L(NS_SYM_L, 1, ndim, 1, 5);

FArray1<double> div1(NS_SYM(div1), 1, ndim);
FArray1<double> p_s(NS_SYM(p_s), 1, ndim);
FArray1<double> r_s(NS_SYM(r_s), 1, ndim);
FArray1<double> r2_s(NS_SYM(r2_s), 1, ndim);
FArray1<double> q_s(NS_SYM(q_s), 1, ndim);
FArray1<double> s_s(NS_SYM(s_s), 1, ndim);
FArray1<double> x1(NS_SYM(x1), 1, ndim);

FArray1<double> xx(NS_SYM(xx), 1, nx);
FArray1<double> yy(NS_SYM(yy), 1, ny);

}  // namespace ns
