// C++ port of calcul_pressure_field.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void calcul_pressure_field() {
    int k = 1;
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            pre(i, j) = x1(k);
            k = k + 1;
        }
    }

    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            p(i, j) = pre(i, j);
        }
    }
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces calcul_pressure_field.f90 inside the Fortran program
extern "C" void calcul_pressure_field_() { ns::calcul_pressure_field(); }
#endif
