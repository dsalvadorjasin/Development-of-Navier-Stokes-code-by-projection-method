// C++ port of discretisation_upwind.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void discretisation_intermediary_velocity_upwind() {
    // u_star calculation (x component)
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            const double u_tilde_x1 = 0.5 * (u(i, j) + u(i + 1, j));
            const double u_tilde_x2 = 0.5 * (u(i - 1, j) + u(i, j));
            const double v_tilde_x1 = 0.5 * (v(i, j) + v(i + 1, j));
            const double v_tilde_x2 = 0.5 * (v(i, j - 1) + v(i + 1, j - 1));

            double a, b, c, d;

            if (u_tilde_x1 >= 0) {
                a = u(i, j);
            } else {
                a = u(i + 1, j);
            }

            if (u_tilde_x2 >= 0) {
                b = u(i - 1, j);
            } else {
                b = u(i, j);
            }

            if (v_tilde_x1 >= 0) {
                c = u(i, j);
            } else {
                c = u(i, j + 1);
            }

            if (v_tilde_x2 >= 0) {
                d = u(i, j - 1);
            } else {
                d = u(i, j);
            }

            u_star(i, j) = u(i, j) - dt * ddx * (u_tilde_x1 * a - u_tilde_x2 * b) -
                           dt * ddy * (v_tilde_x1 * c - v_tilde_x2 * d) +
                           nu * dt * (ddx * ddx) * (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j)) +
                           nu * dt * (ddy * ddy) * (u(i, j + 1) - 2 * u(i, j) + u(i, j - 1));
        }
    }

    // v_star calculation (y component)
    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx; ++i) {
            const double u_tilde_y1 = 0.5 * (u(i, j) + u(i, j + 1));
            const double u_tilde_y2 = 0.5 * (u(i - 1, j) + u(i - 1, j + 1));
            const double v_tilde_y1 = 0.5 * (v(i, j) + v(i, j + 1));
            const double v_tilde_y2 = 0.5 * (v(i, j - 1) + v(i, j));

            double e, f, g, h;

            if (u_tilde_y1 >= 0) {
                e = v(i, j);
            } else {
                e = v(i + 1, j);
            }

            if (u_tilde_y2 >= 0) {
                f = v(i - 1, j);
            } else {
                f = v(i, j);
            }

            if (v_tilde_y1 >= 0) {
                g = v(i, j);
            } else {
                g = v(i, j + 1);
            }

            if (v_tilde_y2 >= 0) {
                h = v(i, j - 1);
            } else {
                h = v(i, j);
            }

            v_star(i, j) = v(i, j) - dt * ddx * (u_tilde_y1 * e - u_tilde_y2 * f) -
                           dt * ddy * (v_tilde_y1 * g - v_tilde_y2 * h) +
                           nu * dt * (ddx * ddx) * (v(i + 1, j) - 2 * v(i, j) + v(i - 1, j)) +
                           nu * dt * (ddy * ddy) * (v(i, j + 1) - 2 * v(i, j) + v(i, j - 1));
        }
    }
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces discretisation_upwind.f90 inside the Fortran program
extern "C" void discretisation_intermediary_velocity_upwind_() { ns::discretisation_intermediary_velocity_upwind(); }
#endif
