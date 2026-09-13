// C++ port of discretisation_centre.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void discretisation_intermediary_velocity_centre() {
    double u_tilde_x1, u_tilde_x2, v_tilde_x1, v_tilde_x2;
    double u_tilde_y1, u_tilde_y2, v_tilde_y1, v_tilde_y2;

    // u_star calculation (x component)
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            u_tilde_x1 = 0.5 * (u(i, j) + u(i + 1, j));
            u_tilde_x2 = 0.5 * (u(i - 1, j) + u(i, j));
            v_tilde_x1 = 0.5 * (v(i, j) + v(i + 1, j));
            v_tilde_x2 = 0.5 * (v(i, j - 1) + v(i + 1, j - 1));

            u_star(i, j) = u(i, j)
                - dt * ddx * 0.5 * (u_tilde_x1 * (u(i, j) + u(i + 1, j)) - u_tilde_x2 * (u(i, j) + u(i - 1, j)))
                - dt * ddy * 0.5 * (v_tilde_x1 * (u(i, j) + u(i, j + 1)) - v_tilde_x2 * (u(i, j) + u(i, j - 1)))
                + nu * dt * (ddx * ddx) * (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j))
                + nu * dt * (ddy * ddy) * (u(i, j + 1) - 2 * u(i, j) + u(i, j - 1));
        }
    }

    // v_star calculation (y component)
    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx; ++i) {
            u_tilde_y1 = 0.5 * (u(i, j) + u(i, j + 1));
            u_tilde_y2 = 0.5 * (u(i - 1, j) + u(i - 1, j + 1));
            v_tilde_y1 = 0.5 * (v(i, j) + v(i, j + 1));
            v_tilde_y2 = 0.5 * (v(i, j - 1) + v(i, j));

            v_star(i, j) = v(i, j)
                - dt * ddx * 0.5 * (u_tilde_y1 * (v(i, j) + v(i + 1, j)) - u_tilde_y2 * (v(i, j) + v(i - 1, j)))
                - dt * ddy * 0.5 * (v_tilde_y1 * (v(i, j) + v(i, j + 1)) - v_tilde_y2 * (v(i, j) + v(i, j - 1)))
                + nu * dt * (ddx * ddx) * (v(i + 1, j) - 2 * v(i, j) + v(i - 1, j))
                + nu * dt * (ddy * ddy) * (v(i, j + 1) - 2 * v(i, j) + v(i, j - 1));
        }
    }
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces discretisation_centre.f90 inside the Fortran program
extern "C" void discretisation_intermediary_velocity_centre_() { ns::discretisation_intermediary_velocity_centre(); }
#endif
