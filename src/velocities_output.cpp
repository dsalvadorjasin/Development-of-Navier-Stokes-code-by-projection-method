// C++ port of velocities_output.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void velocities_output() { not_ported("velocities_output"); }

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces velocities_output.f90 inside the Fortran program
extern "C" void velocities_output_() { ns::velocities_output(); }
#endif
