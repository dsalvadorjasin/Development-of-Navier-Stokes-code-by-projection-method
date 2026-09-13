// C++ port of ensight_cavity.f90: Ensight Gold binary output
// (write_result_ensight, WriteRectEnsightGeo, WriteEnsightVar, EnsightCase).
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void write_result_ensight(FArray1<double>& xEnsight, FArray1<double>& yEnsight,
                          FArray2<double>& u, FArray2<double>& v, FArray2<double>& omega,
                          FArray2<double>& press, int nx, int ny, int nz, int istep, int isto,
                          int nstep) {
    (void)xEnsight; (void)yEnsight; (void)u; (void)v; (void)omega; (void)press; (void)nx;
    (void)ny; (void)nz; (void)istep; (void)isto; (void)nstep;
    not_ported("write_result_ensight");
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces ensight_cavity.f90 inside the Fortran program
extern "C" void write_result_ensight_(double* x, double* y, double* u, double* v, double* omega,
                                      double* press, int* nx, int* ny, int* nz, int* istep,
                                      int* isto, int* nstep) {
    ns::FArray1<double> x_v(x, 1, *nx), y_v(y, 1, *ny);
    ns::FArray2<double> u_v(u, 1, *nx, 1, *ny), v_v(v, 1, *nx, 1, *ny),
        omega_v(omega, 1, *nx, 1, *ny), press_v(press, 1, *nx, 1, *ny);
    ns::write_result_ensight(x_v, y_v, u_v, v_v, omega_v, press_v, *nx, *ny, *nz, *istep, *isto,
                             *nstep);
}
#endif
