// C++ port of mesh_creation.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void mesh_creation() { not_ported("mesh_creation"); }

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces mesh_creation.f90 inside the Fortran program
extern "C" void mesh_creation_() { ns::mesh_creation(); }
#endif
