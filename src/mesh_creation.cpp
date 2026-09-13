// C++ port of mesh_creation.f90
#include "routines.hpp"
#include "variables.hpp"

#include <fstream>
#include <iomanip>

namespace ns {

void mesh_creation() {
    for (int i = 1; i <= nx; ++i) {
        xx(i) = (i - 0.5) * dx;
    }

    for (int j = 1; j <= ny; ++j) {
        yy(j) = (j - 0.5) * dy;
    }

    std::ofstream output("coordinates.txt", std::ios::trunc);
    output << std::setw(6) << nx << std::setw(6) << ny << '\n';
    output << std::fixed << std::setprecision(3);
    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            output << std::setw(6) << xx(i) << "      "
                   << std::setw(6) << yy(j) << '\n';
        }
    }
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces mesh_creation.f90 inside the Fortran program
extern "C" void mesh_creation_() { ns::mesh_creation(); }
#endif
