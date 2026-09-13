// C++ port of calcul_new_velocity.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void calcul_new_velocity() { not_ported("calcul_new_velocity"); }

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces calcul_new_velocity.f90 inside the Fortran program
extern "C" void calcul_new_velocity_() { ns::calcul_new_velocity(); }
#endif
