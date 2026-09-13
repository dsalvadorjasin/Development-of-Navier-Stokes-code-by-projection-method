// C++ port of matgen_cavity.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void matgen_cavity(FArray2<double>& coef, FArray1<int>& jcoef, int nx, int ny, int ndim,
                   int mdim, double dx, double dy) {
    (void)coef; (void)jcoef; (void)nx; (void)ny; (void)ndim; (void)mdim; (void)dx; (void)dy;
    not_ported("matgen_cavity");
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces matgen_cavity.f90 inside the Fortran program
extern "C" void matgen_cavity_(double* coef, int* jcoef, int* nx, int* ny, int* ndim, int* mdim,
                               double* dx, double* dy) {
    ns::FArray2<double> coef_v(coef, 1, *ndim, 1, *mdim);
    ns::FArray1<int> jcoef_v(jcoef, 1, *mdim);
    ns::matgen_cavity(coef_v, jcoef_v, *nx, *ny, *ndim, *mdim, *dx, *dy);
}
#endif
