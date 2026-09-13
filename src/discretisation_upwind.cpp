// C++ port of discretisation_upwind.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void discretisation_intermediary_velocity_upwind() { not_ported("discretisation_intermediary_velocity_upwind"); }

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces discretisation_upwind.f90 inside the Fortran program
extern "C" void discretisation_intermediary_velocity_upwind_() { ns::discretisation_intermediary_velocity_upwind(); }
#endif
