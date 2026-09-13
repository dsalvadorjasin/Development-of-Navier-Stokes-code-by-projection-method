// C++ port of initial_conditions.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void initial_conditions() { not_ported("initial_conditions"); }

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces initial_conditions.f90 inside the Fortran program
extern "C" void initial_conditions_() { ns::initial_conditions(); }
#endif
