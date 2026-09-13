// C++ port of velocities_output.f90
#include <cstdio>

#include "routines.hpp"
#include "variables.hpp"

namespace ns {

namespace {

// Fortran `write(unit,'(F1.5)') x`: no value fits in a 1-character field, so
// gfortran prints a single `*` per record whatever the value is.
void write_f1_5(std::FILE* f, double /*x*/) { std::fputs("*\n", f); }

}  // namespace

void velocities_output() {
    std::FILE* f9 = std::fopen("u_velocity.txt", "w");
    std::fputs(" u_vertical_line\n", f9);  // list-directed write: leading blank
    for (int j = 1; j <= ny; ++j) write_f1_5(f9, u((nx - 1) / 2, j));
    std::fclose(f9);

    std::FILE* f10 = std::fopen("v_velocity.txt", "w");
    std::fputs(" v_horizontal_line\n", f10);
    for (int i = 1; i <= nx; ++i) write_f1_5(f10, v(i, (ny - 1) / 2));
    std::fclose(f10);

    std::FILE* f11 = std::fopen("omega_boundary.txt", "w");
    std::fputs(" omega_boundary\n", f11);
    for (int i = 2; i <= nx; ++i) write_f1_5(f11, rot(i, ny));
    std::fclose(f11);
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces velocities_output.f90 inside the Fortran program
extern "C" void velocities_output_() { ns::velocities_output(); }
#endif
