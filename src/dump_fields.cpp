// C++ port of dump_fields.f90: optional full-precision dump of the raw field
// arrays (env NS_DUMP_FIELDS) used to validate the port against the Fortran
// reference (see tools/compare_reference.py).
#include <cstdio>
#include <cstdlib>

#include "routines.hpp"
#include "variables.hpp"

namespace ns {

namespace {

void dump_array(std::FILE* f, const char* name, const FArray2<double>& a) {
    std::fprintf(f, "%s %5d %5d %5d %5d\n", name, a.lo1(), a.hi1(), a.lo2(), a.hi2());
    for (int j = a.lo2(); j <= a.hi2(); ++j)
        for (int i = a.lo1(); i <= a.hi1(); ++i) std::fprintf(f, "%25.17E\n", a(i, j));
}

}  // namespace

void dump_fields() {
    if (std::getenv("NS_DUMP_FIELDS") == nullptr) return;
    if (!(istep <= 3 || istep % isto == 0)) return;

    char fname[32];
    std::snprintf(fname, sizeof fname, "fields_%06d.txt", istep);
    std::FILE* f = std::fopen(fname, "w");
    if (!f) return;
    std::fprintf(f, "istep %8d\n", istep);
    std::fprintf(f, "time %25.17E\n", time);
    std::fprintf(f, "dt %25.17E\n", dt);
    dump_array(f, "u", u);
    dump_array(f, "v", v);
    dump_array(f, "p", p);
    dump_array(f, "u_star", u_star);
    dump_array(f, "v_star", v_star);
    dump_array(f, "rot", rot);
    dump_array(f, "div", div);
    std::fclose(f);
}

}  // namespace ns

#ifdef NS_HYBRID
extern "C" void dump_fields_() { ns::dump_fields(); }
#endif
