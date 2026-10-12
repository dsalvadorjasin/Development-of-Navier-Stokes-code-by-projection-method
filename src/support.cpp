#include <cstdio>
#include <cstdlib>

#include "routines.hpp"

namespace ns {

void io_fail(const char* filename) {
    std::fprintf(stderr, "error: cannot write output file '%s'\n", filename);
    std::exit(2);
}

void not_ported(const char* routine) {
    std::fprintf(stderr, "error: routine '%s' has not been ported to C++ yet\n", routine);
    std::abort();
}

}  // namespace ns
