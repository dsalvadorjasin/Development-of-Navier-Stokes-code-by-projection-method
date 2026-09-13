// C++ port of timestep.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void timestep() { not_ported("timestep"); }

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces timestep.f90 inside the Fortran program
extern "C" void timestep_() { ns::timestep(); }
#endif
