// C++ port of discretisation_centre.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void discretisation_intermediary_velocity_centre() { not_ported("discretisation_intermediary_velocity_centre"); }

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces discretisation_centre.f90 inside the Fortran program
extern "C" void discretisation_intermediary_velocity_centre_() { ns::discretisation_intermediary_velocity_centre(); }
#endif
