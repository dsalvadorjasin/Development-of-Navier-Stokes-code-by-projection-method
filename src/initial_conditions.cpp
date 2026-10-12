// C++ port of initial_conditions.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void initial_conditions() {
    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i) {
            u(i, j) = 0;
            v(i, j) = 0;
            u1(i, j) = 0;
            v1(i, j) = 0;
            p(i, j) = 0;
            u_star(i, j) = 0;
            v_star(i, j) = 0;
        }
    }
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces initial_conditions.f90 inside the Fortran program
extern "C" void initial_conditions_() { ns::initial_conditions(); }
#endif
