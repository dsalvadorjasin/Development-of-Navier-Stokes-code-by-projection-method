// C++ port of calcul_pressure_field.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void calcul_pressure_field() { not_ported("calcul_pressure_field"); }

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces calcul_pressure_field.f90 inside the Fortran program
extern "C" void calcul_pressure_field_() { ns::calcul_pressure_field(); }
#endif
