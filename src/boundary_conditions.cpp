// C++ port of boundary_conditions.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void boundary_conditions() {
    // Velocity boundary conditions

    // Left vertical wall
    for (int j = 0; j <= ny + 1; ++j) {
        u(0, j) = 0;
        v(0, j) = -v(1, j);
        u1(0, j) = 0;
        v1(0, j) = -v1(1, j);
        u_star(0, j) = 0;
        v_star(0, j) = -v_star(1, j);
    }

    // Right vertical wall
    for (int j = 0; j <= ny + 1; ++j) {
        u(nx, j) = 0;
        v(nx + 1, j) = -v(nx, j);
        u1(nx, j) = 0;
        v1(nx + 1, j) = -v1(nx, j);
        u_star(nx, j) = 0;
        v_star(nx + 1, j) = -v_star(nx, j);
    }

    // Bottom horizontal wall
    for (int i = 0; i <= nx + 1; ++i) {
        u(i, 0) = -u(i, 1);
        v(i, 0) = 0;
        u1(i, 0) = -u1(i, 1);
        v1(i, 0) = 0;
        u_star(i, 0) = -u_star(i, 1);
        v_star(i, 0) = 0;
    }

    // Top horizontal wall
    for (int i = 0; i <= nx + 1; ++i) {
        u(i, ny + 1) = 2 - u(i, ny);
        v(i, ny) = 0;
        u1(i, ny + 1) = 2 - u1(i, ny);
        v1(i, ny) = 0;
        u_star(i, ny + 1) = 2 - u_star(i, ny);
        v_star(i, ny) = 0;
    }

    // Pressure boundary conditions

    // Left vertical wall
    for (int j = 0; j <= ny + 1; ++j) {
        p(0, j) = p(1, j);
    }

    // Right vertical wall
    for (int j = 0; j <= ny + 1; ++j) {
        p(nx + 1, j) = p(nx, j);
    }

    // Bottom horizontal wall
    for (int i = 0; i <= nx + 1; ++i) {
        p(i, 0) = p(i, 1);
    }

    // Top horizontal wall
    for (int i = 0; i <= nx + 1; ++i) {
        p(i, ny + 1) = p(i, ny);
    }
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces boundary_conditions.f90 inside the Fortran program
extern "C" void boundary_conditions_() { ns::boundary_conditions(); }
#endif
