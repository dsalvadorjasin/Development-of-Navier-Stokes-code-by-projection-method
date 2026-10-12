// C++ port of reading_data.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void reading_data() {
    dx = 1. / double(nx);
    dy = 1. / double(ny);
    ddx = 1. / dx;
    ddy = 1. / dy;
    nu1 = 1. / nu;
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces reading_data.f90 inside the Fortran program
extern "C" void reading_data_() { ns::reading_data(); }
#endif
