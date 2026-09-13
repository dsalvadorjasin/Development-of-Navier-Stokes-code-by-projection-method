// C++ port of calcul_2nd_membre.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void calcul_2nd_membre() { not_ported("calcul_2nd_membre"); }

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces calcul_2nd_membre.f90 inside the Fortran program
extern "C" void calcul_2nd_membre_() { ns::calcul_2nd_membre(); }
#endif
