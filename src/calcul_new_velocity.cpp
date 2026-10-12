// C++ port of calcul_new_velocity.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void calcul_new_velocity() {
    int i, j;

    // x direction
    for (j = 1; j <= ny; ++j) {
        for (i = 1; i <= nx - 1; ++i) {
            u1(i, j) = u_star(i, j) - dt * ddx * (p(i + 1, j) - p(i, j));
        }
    }

    // y direction
    for (j = 1; j <= ny - 1; ++j) {
        for (i = 1; i <= nx; ++i) {
            v1(i, j) = v_star(i, j) - dt * ddy * (p(i, j + 1) - p(i, j));
        }
    }

    // centring to (1:nx,1:ny)
    for (j = 1; j <= ny; ++j) {
        for (i = 1; i <= nx; ++i) {
            uc(i, j) = 0.5 * (u1(i, j) + u1(i - 1, j));
            vc(i, j) = 0.5 * (v1(i, j) + v1(i, j - 1));
        }
    }

    // rotational: corners
    w(1, 1) = 0;
    w(nx, 1) = 0;
    w(1, ny) = 0;
    w(nx, ny) = 0;

    // general calculation; leaves j == ny+1, as in the Fortran
    for (j = 2; j <= ny; ++j) {
        for (i = 2; i <= nx; ++i) {
            w(i, j) = ddx * (v1(i + 1, j) - v1(i, j)) - ddy * (u1(i, j + 1) - u1(i, j));
        }
    }

    // top and bottom walls (u1(i,j) here reads u1(i,ny+1), Fortran quirk kept)
    for (i = 2; i <= nx - 1; ++i) {
        w(i, 1) = ddx * (v1(i + 1, 1) - v1(i, 1)) - ddy * (u1(i, 2) - u1(i, 1));
        w(i, ny) = ddx * (v1(i + 1, ny) - v1(i, ny)) - ddy * (u1(i, ny + 1) - u1(i, j));
    }

    // left and right walls
    for (j = 2; j <= ny - 1; ++j) {
        w(1, j) = ddx * (v1(2, j) - v1(1, j)) - ddy * (u1(1, j + 1) - u1(1, j));
        w(nx, j) = ddx * (v1(nx + 1, j) - v1(nx, j)) - ddy * (u1(nx, j + 1) - u1(nx, j));
    }

    // centred rotational
    for (j = 2; j <= ny; ++j) {
        for (i = 2; i <= nx; ++i) {
            rot(i, j) = 0.25 * (w(i, j) + w(i - 1, j) + w(i, j - 1) + w(i - 1, j - 1));
        }
    }
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces calcul_new_velocity.f90 inside the Fortran program
extern "C" void calcul_new_velocity_() { ns::calcul_new_velocity(); }
#endif
