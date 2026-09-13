#include <cstdio>
#include <cstdlib>

#include "routines.hpp"

namespace ns {

void not_ported(const char* routine) {
    std::fprintf(stderr, "error: routine '%s' has not been ported to C++ yet\n", routine);
    std::abort();
}

}  // namespace ns
