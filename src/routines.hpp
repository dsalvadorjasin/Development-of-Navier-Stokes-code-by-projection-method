// Declarations of every routine of the solver (one .cpp per original .f90).
// Routines that take no arguments operate on the globals of variables.hpp,
// exactly like their Fortran counterparts operate on `use variables`.
#pragma once

#include "fortran_array.hpp"

namespace ns {

// setup                                                 (original file)
void reading_data();                                  // reading_data.f90
void mesh_creation();                                 // mesh_creation.f90
void initial_conditions();                            // initial_conditions.f90
void boundary_conditions();                           // boundary_conditions.f90
void timestep();                                      // timestep.f90

// momentum discretisation
void discretisation_intermediary_velocity_upwind();   // discretisation_upwind.f90
void discretisation_intermediary_velocity_centre();   // discretisation_centre.f90

// pressure / projection step
void calcul_2nd_membre();                             // calcul_2nd_membre.f90
void matgen_cavity(FArray2<double>& coef, FArray1<int>& jcoef, int nx, int ny, int ndim,
                   int mdim, double dx, double dy);   // matgen_cavity.f90
void ICCG2(FArray2<double>& coef, FArray1<int>& jcoef, FArray2<double>& L, FArray1<int>& Ldiag,
           FArray1<double>& b, FArray1<double>& x, int n, int m, double eps, FArray1<double>& p,
           FArray1<double>& r, FArray1<double>& r2, FArray1<double>& q, FArray1<double>& s,
           int itmax);                                // solver.f90
void calcul_pressure_field();                         // calcul_pressure_field.f90
void calcul_new_velocity();                           // calcul_new_velocity.f90

// output
void velocities_output();                             // velocities_output.f90
void write_result_ensight(FArray1<double>& xEnsight, FArray1<double>& yEnsight,
                          FArray2<double>& u, FArray2<double>& v, FArray2<double>& omega,
                          FArray2<double>& press, int nx, int ny, int nz, int istep, int isto,
                          int nstep);                 // ensight_cavity.f90
void dump_fields();                                   // dump_fields.f90 (validation aid)

// Reports a failed open/write/close of an output file and exits with status 2,
// like a gfortran I/O runtime error (the Fortran has no iostat= checks).
[[noreturn]] void io_fail(const char* filename);

// Called by the stubs of routines that have not been ported yet.
[[noreturn]] void not_ported(const char* routine);

}  // namespace ns
