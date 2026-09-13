// Lightweight Fortran-style array views used throughout the C++ port.
//
// They keep the ORIGINAL Fortran index ranges (e.g. u(0:nx+1,0:ny+1),
// coef(1:ndim,1:mdim)) and the Fortran COLUMN-MAJOR layout (first index
// varies fastest), so every Fortran expression `a(i,j)` becomes `a(i,j)` in
// C++ without any index arithmetic.  This is also what makes the arrays
// interchangeable with the gfortran module storage in the hybrid build.
//
// Compile with -DNS_BOUNDS_CHECK to assert on out-of-range access.
#pragma once

#include <cstdio>
#include <cstdlib>

#ifdef NS_BOUNDS_CHECK
#define NS_CHECK(cond, what)                                                        \
    do {                                                                            \
        if (!(cond)) {                                                              \
            std::fprintf(stderr, "array bounds error: %s (%s:%d)\n", what, __FILE__, \
                         __LINE__);                                                 \
            std::abort();                                                           \
        }                                                                           \
    } while (0)
#else
#define NS_CHECK(cond, what) ((void)0)
#endif

namespace ns {

template <class T>
class FArray1 {
public:
    FArray1(T* data, int lo, int hi) : data_(data), lo_(lo), hi_(hi) {}

    T& operator()(int i) const {
        NS_CHECK(i >= lo_ && i <= hi_, "FArray1 index");
        return data_[i - lo_];
    }
    int lo() const { return lo_; }
    int hi() const { return hi_; }
    int size() const { return hi_ - lo_ + 1; }
    T* data() const { return data_; }

    void fill(T v) const {
        for (int n = 0; n < size(); ++n) data_[n] = v;
    }
    // Fortran maxval()
    T maxval() const {
        T m = data_[0];
        for (int n = 1; n < size(); ++n)
            if (data_[n] > m) m = data_[n];
        return m;
    }

private:
    T* data_;
    int lo_, hi_;
};

template <class T>
class FArray2 {
public:
    FArray2(T* data, int lo1, int hi1, int lo2, int hi2)
        : data_(data), lo1_(lo1), hi1_(hi1), lo2_(lo2), hi2_(hi2) {}

    // column-major: a(i,j) -> data[(i-lo1) + (j-lo2)*n1]
    T& operator()(int i, int j) const {
        NS_CHECK(i >= lo1_ && i <= hi1_, "FArray2 first index");
        NS_CHECK(j >= lo2_ && j <= hi2_, "FArray2 second index");
        return data_[(i - lo1_) + (j - lo2_) * n1()];
    }
    int lo1() const { return lo1_; }
    int hi1() const { return hi1_; }
    int lo2() const { return lo2_; }
    int hi2() const { return hi2_; }
    int n1() const { return hi1_ - lo1_ + 1; }
    int n2() const { return hi2_ - lo2_ + 1; }
    int size() const { return n1() * n2(); }
    T* data() const { return data_; }

    void fill(T v) const {
        for (int n = 0; n < size(); ++n) data_[n] = v;
    }
    // Fortran maxval()
    T maxval() const {
        T m = data_[0];
        for (int n = 1; n < size(); ++n)
            if (data_[n] > m) m = data_[n];
        return m;
    }

private:
    T* data_;
    int lo1_, hi1_, lo2_, hi2_;
};

}  // namespace ns
