// C++ port of solver.f90: ICCG2 (conjugate gradient + incomplete LU fill-in 2
// preconditioner) and its helpers resLU_ell5, matmul_ell, scal, norm1, norm2.
//
// Storage: coef(1:n,1:m) holds the symmetric matrix in ELL/diagonal form,
// coef(i,j) being the entry (i, i+jcoef(j)) (jcoef(1)=0 is the diagonal).
// L(1:n,1:5) holds the incomplete factor with the five offsets Ldiag(1:5) =
// (0, jcoef(2), jcoef(3), jcoef(3)-1, jcoef(3)-2); L(i,j) is the entry
// (i, i+Ldiag(j)).  All arrays are FArray views with the Fortran bounds, so
// every index expression below is a literal transcription of the Fortran.
#include <cmath>
#include <cstdio>

#include "routines.hpp"
#include "variables.hpp"

namespace ns {

// x = L^-T L^-1 y, forward then backward sweep over the five diagonals of L
void resLU_ell5(FArray1<double>& x, FArray2<double>& L, FArray1<int>& Ldiag, FArray1<double>& y,
                int n) {
    int i, col1, col2, col3, col4;

    col1 = Ldiag(2);
    col2 = Ldiag(5);
    col3 = Ldiag(4);
    col4 = Ldiag(3);

    for (i = 1; i <= col1; ++i) {
        x(i) = y(i);
    }

    for (i = col1 + 1; i <= col2; ++i) {
        x(i) = y(i) - L(i - col1, 2) * x(i - col1);
    }

    for (i = col2 + 1; i <= col3; ++i) {
        x(i) = y(i) - L(i - col1, 2) * x(i - col1) - L(i - col2, 5) * x(i - col2);
    }

    for (i = col3 + 1; i <= col4; ++i) {
        x(i) = y(i) - L(i - col1, 2) * x(i - col1) - L(i - col2, 5) * x(i - col2)
               - L(i - col3, 4) * x(i - col3);
    }

    for (i = col4 + 1; i <= n; ++i) {
        x(i) = y(i) - L(i - col1, 2) * x(i - col1) - L(i - col2, 5) * x(i - col2)
               - L(i - col3, 4) * x(i - col3) - L(i - col4, 3) * x(i - col4);
    }

    for (i = n - col1; i >= n - col2 + 1; --i) {
        x(i) = x(i) - L(i, 2) * x(i + col1);
    }

    for (i = n - col2; i >= n - col3 + 1; --i) {
        x(i) = x(i) - L(i, 2) * x(i + col1) - L(i, 5) * x(i + col2);
    }

    for (i = n - col3; i >= n - col4 + 1; --i) {
        x(i) = x(i) - L(i, 2) * x(i + col1) - L(i, 5) * x(i + col2) - L(i, 4) * x(i + col3);
    }

    for (i = n - col4; i >= 1; --i) {
        x(i) = x(i) - L(i, 2) * x(i + col1) - L(i, 5) * x(i + col2) - L(i, 4) * x(i + col3)
               - L(i, 3) * x(i + col4);
    }
}

// x = A y for the symmetric matrix A stored by diagonals (coef, jcoef)
void matmul_ell(FArray1<double>& x, FArray2<double>& coef, FArray1<int>& jcoef, FArray1<double>& y,
                int n, int m) {
    int i, j, col;

    for (i = 1; i <= n; ++i) {
        x(i) = coef(i, 1) * y(i);
    }

    for (j = 2; j <= m; ++j) {
        col = jcoef(j);
        for (i = 1; i <= n - col; ++i) {
            x(i) = x(i) + coef(i, j) * y(i + col);
            x(i + col) = x(i + col) + coef(i, j) * y(i);
        }
    }
}

double scal(FArray1<double>& x, FArray1<double>& y, int n) {
    int i;
    double res;

    res = 0.;
    for (i = 1; i <= n; ++i) res = res + x(i) * y(i);
    return res;
}

double norm1(FArray1<double>& x, int n) {
    int i;
    double res;

    res = 0.;
    for (i = 1; i <= n; ++i) res = res + x(i) * x(i);
    return res;
}

double norm2(FArray1<double>& x, int n) {
    int i;
    double res;

    res = 0.;
    for (i = 1; i <= n; ++i) res = res + x(i) * x(i);
    res = std::sqrt(res);
    return res;
}

// ATTENTION: coef, jcoef and b are modified in this routine (as in the Fortran)
void ICCG2(FArray2<double>& coef, FArray1<int>& jcoef, FArray2<double>& L, FArray1<int>& Ldiag,
           FArray1<double>& b, FArray1<double>& x, int n, int m, double eps, FArray1<double>& p,
           FArray1<double>& r, FArray1<double>& r2, FArray1<double>& q, FArray1<double>& s,
           int itmax) {
    double alpha, beta, nu, mu, norm0, norm;
    int i, j, col;

    std::printf(" itmax= %11d\n", itmax);
    std::fflush(stdout);

    //-----------------------------------------------------------!
    //   PRECONDITIONING MATRIX CALCULUS + SCALING               !
    //-----------------------------------------------------------!
    for (i = 1; i <= 3; ++i) Ldiag(i) = jcoef(i);
    Ldiag(4) = Ldiag(3) - 1;
    Ldiag(5) = Ldiag(3) - 2;

    L.fill(0.);

    for (i = 1; i <= n; ++i) {
        L(i, 1) = coef(i, 1);
        for (j = 2; j <= 5; ++j) {
            col = Ldiag(j);
            if (col < i) L(i, 1) = L(i, 1) - L(i - col, j) * L(i - col, j) * s(i - col);
        }
        s(i) = 1. / L(i, 1);  // Temporary variable to avoid divisions by L(i,1)

        L(i, 2) = coef(i, 2);
        for (j = 4; j <= 5; ++j) {
            col = Ldiag(j);
            if (col < i) L(i, 2) = L(i, 2) - L(i - col, j) * L(i - col, j - 1) * s(i - col);
        }

        L(i, 3) = coef(i, 3);

        if (i > 1) {
            L(i, 4) = -L(i - 1, 3) * L(i - 1, 2) * s(i - 1);
            L(i, 5) = -L(i - 1, 4) * L(i - 1, 2) * s(i - 1);
        }
    }

    for (i = 1; i <= n; ++i) L(i, 1) = std::sqrt(s(i));

    for (j = 1; j <= 3; ++j) {
        col = Ldiag(j);
        for (i = 1; i <= n - col; ++i) {
            coef(i, j) = coef(i, j) * L(i, 1) * L(i + col, 1);
        }
    }

    for (i = 1; i <= n; ++i) {
        L(i, 2) = coef(i, 2);
        for (j = 4; j <= 5; ++j) {
            col = Ldiag(j);
            if (col < i) L(i, 2) = L(i, 2) - L(i - col, j) * L(i - col, j - 1);
        }

        L(i, 3) = coef(i, 3);

        if (i > 1) {
            L(i, 4) = -L(i - 1, 3) * L(i - 1, 2);
            L(i, 5) = -L(i - 1, 4) * L(i - 1, 2);
        }
    }

    norm0 = norm2(b, n);

    for (i = 1; i <= n; ++i) s(i) = 1. / s(i);

    for (i = 1; i <= n; ++i) b(i) = b(i) * L(i, 1);

    //---------------------------------------------------!
    //    CONJUGUATE GRADIENT                            !
    //---------------------------------------------------!

    matmul_ell(p, coef, jcoef, x, n, m);
    for (i = 1; i <= n; ++i) r(i) = b(i) - p(i);

    resLU_ell5(r2, L, Ldiag, r, n);

    for (i = 1; i <= n; ++i) p(i) = r2(i);
    nu = scal(r, r2, n);

    norm = 0.;
    for (i = 1; i <= n; ++i) norm = norm + r(i) * r(i) * s(i);
    norm = std::sqrt(norm) / norm0;

    j = 0;

    while (norm > eps && j < itmax) {
        j = j + 1;

        matmul_ell(q, coef, jcoef, p, n, m);

        alpha = nu / scal(p, q, n);

        for (i = 1; i <= n; ++i) x(i) = x(i) + alpha * p(i);
        for (i = 1; i <= n; ++i) r(i) = r(i) - alpha * q(i);

        resLU_ell5(r2, L, Ldiag, r, n);

        mu = scal(r, r2, n);

        beta = mu / nu;

        for (i = 1; i <= n; ++i) p(i) = r2(i) + beta * p(i);

        nu = mu;

        norm = 0.;
        for (i = 1; i <= n; ++i) norm = norm + r(i) * r(i) * s(i);
        norm = std::sqrt(norm) / norm0;
    }

    //---------------------------------------------------!
    //    SOLUTION SCALING                               !
    //---------------------------------------------------!

    for (i = 1; i <= n; ++i) x(i) = x(i) * L(i, 1);

    // the iteration count must be the LAST token of the line (see PORTING.md);
    // flushed line by line so it cannot interleave with libgfortran's stdout buffer
    if (j >= itmax) {
        std::printf(" non convergence apres = norm=%24.17G it=%11d\n", norm, j);
    } else {
        std::printf("  Nombre Iterations ICCG2 ( Fill-In 2 ) = %11d\n", j);
    }
    std::fflush(stdout);
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces solver.f90 inside the Fortran program
extern "C" void iccg2_(double* coef, int* jcoef, double* L, int* Ldiag, double* b, double* x,
                       int* n, int* m, double* eps, double* p, double* r, double* r2, double* q,
                       double* s, int* itmax) {
    ns::FArray2<double> coef_v(coef, 1, *n, 1, *m);
    ns::FArray1<int> jcoef_v(jcoef, 1, *m);
    ns::FArray2<double> L_v(L, 1, *n, 1, 5);
    ns::FArray1<int> Ldiag_v(Ldiag, 1, 5);
    ns::FArray1<double> b_v(b, 1, *n), x_v(x, 1, *n), p_v(p, 1, *n), r_v(r, 1, *n),
        r2_v(r2, 1, *n), q_v(q, 1, *n), s_v(s, 1, *n);
    ns::ICCG2(coef_v, jcoef_v, L_v, Ldiag_v, b_v, x_v, *n, *m, *eps, p_v, r_v, r2_v, q_v, s_v,
              *itmax);
}
#endif
