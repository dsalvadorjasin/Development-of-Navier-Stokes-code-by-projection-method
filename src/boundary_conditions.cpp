// C++ port of boundary_conditions.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void boundary_conditions() { not_ported("boundary_conditions"); }

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces boundary_conditions.f90 inside the Fortran program
extern "C" void boundary_conditions_() { ns::boundary_conditions(); }
#endif
