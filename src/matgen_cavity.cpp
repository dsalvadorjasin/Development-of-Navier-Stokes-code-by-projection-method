// C++ port of matgen_cavity.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void matgen_cavity(FArray2<double>& coef, FArray1<int>& jcoef, int nx, int ny, int ndim,
                   int mdim, double dx, double dy) {
    (void)ndim; (void)mdim;
    double aa, bb;

    int k = 1;
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            coef(k, 1) = 2. / dx / dx + 2. / dy / dy;
            coef(k, 2) = -1. / dx / dx;
            coef(k, 3) = -1. / dy / dy;
            k = k + 1;
        }
    }

    k = 1;
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            if (i == 1) {
                aa = -1. / dx / dx;
            } else if (i == nx) {
                aa = -1. / dx / dx;
            } else {
                aa = 0.0;
            }

            if (j == 1) {
                bb = -1. / dy / dy;
            } else if (j == ny) {
                bb = -1. / dy / dy;
            } else {
                bb = 0.0;
            }

            coef(k, 1) = coef(k, 1) + aa + bb;
            k = k + 1;
        }
    }

    k = 1;
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            if (i == nx) {
                coef(k, 2) = 0.0;
            }
            if (j == ny) {
                coef(k, 3) = 0.0;
            }
            k = k + 1;
        }
    }

    jcoef(1) = 0;
    jcoef(2) = 1;
    jcoef(3) = nx;
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
