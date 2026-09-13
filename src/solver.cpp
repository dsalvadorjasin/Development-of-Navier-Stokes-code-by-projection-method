// C++ port of solver.f90: ICCG2 (conjugate gradient + incomplete LU fill-in 2
// preconditioner) and its helpers resLU_ell5, matmul_ell, scal, norm1, norm2.
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void ICCG2(FArray2<double>& coef, FArray1<int>& jcoef, FArray2<double>& L, FArray1<int>& Ldiag,
           FArray1<double>& b, FArray1<double>& x, int n, int m, double eps, FArray1<double>& p,
           FArray1<double>& r, FArray1<double>& r2, FArray1<double>& q, FArray1<double>& s,
           int itmax) {
    (void)coef; (void)jcoef; (void)L; (void)Ldiag; (void)b; (void)x; (void)n; (void)m; (void)eps;
    (void)p; (void)r; (void)r2; (void)q; (void)s; (void)itmax;
    not_ported("ICCG2");
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces solver.f90 inside the Fortran program
extern "C" void iccg2_(double* coef, int* jcoef, double* L, int* Ldiag, double* b, double* x,
                       int* n, int* m, double* eps, double* p, double* r, double* r2, double* q,
                       double* s, int* itmax) {
    ns::FArray2<double> coef_v(coef, 1, *n, 1, *m);
    ns::FArray1<int> jcoef_v(jcoef, 1, *m);
    ns::FArray2<double> L_v(L, 1, *n, 1, 5);
    ns::FArray1<int> Ldiag_v(Ldiag, 1, 5);
    ns::FArray1<double> b_v(b, 1, *n), x_v(x, 1, *n), p_v(p, 1, *n), r_v(r, 1, *n),
        r2_v(r2, 1, *n), q_v(q, 1, *n), s_v(s, 1, *n);
    ns::ICCG2(coef_v, jcoef_v, L_v, Ldiag_v, b_v, x_v, *n, *m, *eps, p_v, r_v, r2_v, q_v, s_v,
              *itmax);
}
#endif
