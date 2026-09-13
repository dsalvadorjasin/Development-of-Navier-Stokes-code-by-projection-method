// C++ port of calcul_2nd_membre.f90
#include "routines.hpp"
#include "variables.hpp"

namespace ns {

void calcul_2nd_membre() {
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            div(i, j) = ddt * (ddx * (u_star(i, j) - u_star(i - 1, j)) +
                               ddy * (v_star(i, j) - v_star(i, j - 1)));
        }
    }

    double sum = 0.0;
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            sum = sum + div(i, j);
        }
    }

    sum = sum / double(nx * ny);
    int k = 1;
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            div(i, j) = div(i, j) - sum;
            div1(k) = -div(i, j);
            x1(k) = 0.;
            k = k + 1;
        }
    }
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces calcul_2nd_membre.f90 inside the Fortran program
extern "C" void calcul_2nd_membre_() { ns::calcul_2nd_membre(); }
#endif
