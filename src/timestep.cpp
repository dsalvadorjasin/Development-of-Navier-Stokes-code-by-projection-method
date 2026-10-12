// C++ port of timestep.f90
#include <algorithm>
#include <cmath>

#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void timestep() {
    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i) {
            abs_u(i, j) = std::abs(u(i, j));
            abs_v(i, j) = std::abs(v(i, j));
        }
    }

    max_x = 1. / abs_u.maxval();
    max_y = 1. / abs_v.maxval();

    if (Q == 1) {
        min_u = std::min(dx * max_x, dx * dx * 0.5 * nu1);
        min_v = std::min(dy * max_y, dy * dy * 0.5 * nu1);
    } else if (Q == 2) {
        min_u = std::min(nu * max_x * max_x, dx * dx * 0.5 * nu1);
        min_v = std::min(nu * max_y * max_y, dy * dy * 0.5 * nu1);
    }

    dt = std::min(std::min(min_u, min_v), dt2);
    ddt = 1. / dt;
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces timestep.f90 inside the Fortran program
extern "C" void timestep_() { ns::timestep(); }
#endif
