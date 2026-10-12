// C++ port of velocities_output.f90
#include <cstdio>

#include "routines.hpp"
#include "variables.hpp"

namespace ns {

namespace {

// Fortran `write(unit,'(F1.5)') x`: no value fits in a 1-character field, so
// gfortran prints a single `*` per record whatever the value is.
void write_f1_5(std::FILE* f, double /*x*/) { std::fputs("*\n", f); }

std::FILE* open_output(const char* name) {
    std::FILE* f = std::fopen(name, "w");
    if (!f) io_fail(name);
    return f;
}

void close_output(std::FILE* f, const char* name) {
    const bool failed = std::ferror(f) != 0;
    if (std::fclose(f) != 0 || failed) io_fail(name);
}

}  // namespace

void velocities_output() {
    std::FILE* f9 = open_output("u_velocity.txt");
    std::fputs(" u_vertical_line\n", f9);  // list-directed write: leading blank
    for (int j = 1; j <= ny; ++j) write_f1_5(f9, u((nx - 1) / 2, j));
    close_output(f9, "u_velocity.txt");

    std::FILE* f10 = open_output("v_velocity.txt");
    std::fputs(" v_horizontal_line\n", f10);
    for (int i = 1; i <= nx; ++i) write_f1_5(f10, v(i, (ny - 1) / 2));
    close_output(f10, "v_velocity.txt");

    std::FILE* f11 = open_output("omega_boundary.txt");
    std::fputs(" omega_boundary\n", f11);
    for (int i = 2; i <= nx; ++i) write_f1_5(f11, rot(i, ny));
    close_output(f11, "omega_boundary.txt");
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces velocities_output.f90 inside the Fortran program
extern "C" void velocities_output_() { ns::velocities_output(); }
#endif
