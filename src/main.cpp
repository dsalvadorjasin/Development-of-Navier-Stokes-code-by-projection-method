// C++ port of main.f90: driver / timestep loop of the lid-driven cavity solver.
#include <cstdio>
#include <iostream>

#include "routines.hpp"
#include "variables.hpp"

using namespace ns;

int main() {
    reading_data();         // reading some data for mesh creation and calculation
    mesh_creation();        // creating the mesh for the simulation
    initial_conditions();   // initializing conditions
    boundary_conditions();  // boundary conditions of the domain

    std::cout << " Choice of numerical scheme: 1-Upwind; 2-Centred\n";
    if (!(std::cin >> Q)) {
        std::fprintf(stderr, "error: expected the scheme number (1 or 2) on stdin\n");
        return 1;
    }

    // Parameters for the simulation
    zeta = 1.e-8;   // zeta and itmax for solving pressure matrix in ICCG2
    itmax = 300;
    ns::time = 0.;     // initialize time of simulation
    nstep = 20000;  // number of timesteps for the simulation
    isto = 200;     // data stored every 'isto' steps

    for (istep = 1; istep <= nstep; ++istep) {  // main loop on the timesteps
        if (Q == 1) {  // Upwind scheme
            timestep();
            discretisation_intermediary_velocity_upwind();
        } else if (Q == 2) {  // Centred scheme
            timestep();
            discretisation_intermediary_velocity_centre();
        }

        calcul_2nd_membre();  // Calculate the second member of equation AX=B
        matgen_cavity(coef, jcoef, nx, ny, ndim, mdim, dx, dy);  // Generate coefficients of A
        // CONJUGATE GRADIENT + PRECONDITIONED LU for solving X field (pressure here)
        ICCG2(coef, jcoef, L, Ldiag, div1, x1, ndim, mdim, zeta, p_s, r_s, r2_s, q_s, s_s, itmax);
        calcul_pressure_field();  // update the pressure field
        calcul_new_velocity();    // update velocity field + calculate rotational 'omega'

        // loops to update the main u and v velocity fields
        for (int j = 1; j <= ny; ++j)
            for (int i = 1; i <= nx - 1; ++i) u(i, j) = u1(i, j);

        for (int j = 1; j <= ny - 1; ++j)
            for (int i = 1; i <= nx; ++i) v(i, j) = v1(i, j);

        ns::time = ns::time + dt;  // update in time --> time = time + timestep
        std::printf(" time = %24.17G\n", ns::time);

        boundary_conditions();  // recall boundary conditions to update them

        dump_fields();  // optional full-precision dump (env NS_DUMP_FIELDS)

        if (istep % isto == 0)  // write results in ensight format
            write_result_ensight(xx, yy, uc, vc, rot, pre, nx, ny, nz, istep, isto, nstep);
    }

    velocities_output();  // write velocity and omega on lines at final stage
    return 0;
}
