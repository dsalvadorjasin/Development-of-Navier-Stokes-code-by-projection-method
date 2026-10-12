# Development-of-Navier-Stokes-code-by-projection-method
<p align="justify">This code is a C++ implementation (ported from the original Fortran 90, which is kept in <em>fortran/</em> as the reference) of a 2D flow using projection method with FVM. Navier Stokes equations are solved for velocity and pressure fields. The rotational of velocity is also calculated, and the output data can be visualized with Ensight or Paraview.</p>

## Getting Started

<p align="justify">These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.</p>

### Prerequisites

<p align="justify">The code is written in standard C++17 with no external dependencies; it requires a C++ compiler such as <a href="https://gcc.gnu.org/">g++</a> (or clang++) and GNU make. GFortran is only needed to rebuild the original Fortran program in <em>fortran/</em> (reference oracle, see <a href="PORTING.md">PORTING.md</a>).</p>

In Ubuntu, Mint and Debian you can install them like this:

    sudo apt-get install build-essential          # g++ and make
    sudo apt-get install gfortran                 # optional: reference Fortran build
    
<p align="justify">Alongside the compiler, the open-source platform for data analysis and visualization, <a href="https://www.paraview.org/">ParaView</a>, is also required and installed through the commands:</p>

    sudo apt-get update
    sudo apt-get install paraview

<p align="justify">Paraview is preferred compared to Ensight as it is free and open-source, so that everybody can use it.</p>

For other Linux flavors, OS X and Windows, packages are available at:

https://gcc.gnu.org/install/binaries.html for GCC (g++)
https://www.paraview.org/download/ for ParaView


## File descriptions
<ul>
    <li>'<em>src/*.cpp</em>', '<em>src/*.hpp</em>': the C++ code, one file per routine (<em>main.cpp</em> driver, <em>variables.hpp</em> global state, <em>fortran_array.hpp</em> array views that keep the original index ranges, <em>solver.cpp</em> ICCG2 pressure solver, <em>ensight_cavity.cpp</em> Ensight Gold writer, ...).</li>
    <li>'<em>Makefile</em>' to compile the code (<code>make</code>) and to validate it against the reference (<code>make check</code>).</li>
    <li>'<em>fortran/</em>': the original Fortran 90 program and its Makefile, kept as the reference oracle of the port.</li>
    <li>'<em>reference/</em>': outputs of the Fortran program for both schemes (field dumps, Ensight files, text outputs), '<em>tools/</em>': scripts to regenerate them and to compare a run against them, '<em>PORTING.md</em>': how the port is organised and validated, '<em>docs/REVIEW_solver_ensight.md</em>': review dossier of the solver / Ensight port.</li>
    <li>In the directory '<em>Documents</em>', there are 3 files in PDF format:
        <ul>
            <li>2 papers '<em>cavityflow.pdf</em>' and '<em>ghia82.pdf</em>', which treat of the lid driven cavity problem for incompressible flows.</li>             
            <li>'<em>Program_presentation.pdf</em>' which describes all the steps followed in the code implementation.</li>
        </ul>
    </li>
</ul>

### Running the program

1. Input numerical values for the mesh size and viscosity, in the file 'src/variables.hpp'.
 ```cpp
 // Parameters for reading data to be modified
 constexpr int nx = 11, ny = 11, nz = 1;        // nx and ny correspond to the number of cells in x and y direction
 constexpr double nu = 1e-02, dt2 = 1e-04;      // nu is the fluid viscosity and dt2 is the initial timestep
 ```
2. Simulation parameters to be defined in the main program 'src/main.cpp'.
```cpp
    //!!! Parameters for the simulation
    zeta = 1.e-8;      // zeta and itmax for solving pressure matrix in ICCG2 subroutines
    itmax = 300;
    ns::time = 0.;     // initialize time of simulation
    nstep = 20000;     // number of timesteps for the simulation
    isto = 200;        // data stored every 'isto' steps
```
3. Use the **Makefile** to compile all the files and create the executable (run the command 'make' while being in the main directory of the program).

4. Launch the executable <i>NS_lid_driven_cavity</i>, which will create the mesh and run the discretised calculation on the latter. It asks for the discretisation scheme on standard input: `1` = Upwind, `2` = Centred (e.g. `echo 1 | ./NS_lid_driven_cavity`).

5. Observe the different fields (u, v, p and omega) on the domain using Paraview. Ensight6 and Ensight Gold formats are &emsp;normally supported by Paraview.

6. Remove the created files thanks to the commands 'make clean' and 'make cleanall'.

### Validating against the original Fortran code

The outputs of the original Fortran program (both schemes, default parameters) are stored in `reference/`.
`make check` runs the C++ executable for both schemes and compares every output with them: the field
values (u, v, p, u*, v*, omega, divergence at several timesteps, full double precision), the Ensight
binary files byte for byte and the text outputs exactly. See [PORTING.md](PORTING.md) for the details,
how to regenerate the reference (`tools/make_reference.sh`, needs gfortran) and how to validate a single
routine with the hybrid Fortran/C++ build.
`make check-variants` (needs gfortran) rebuilds both programs with Re=1000 and with a 21x21 grid and
compares the C++ run with the Fortran run for both schemes, with zero tolerance.

### Rendering a video of the velocity field (headless ParaView)

`make_video.py` is a ParaView Python script that loads `EnsightOutput.case`, colours the domain by
velocity magnitude with arrow glyphs, and writes one PNG per stored timestep. The frames are then
assembled into an MP4 with `ffmpeg`.

    sudo apt-get install build-essential paraview ffmpeg   # pvbatch is shipped with paraview

    make                                            # build ./NS_lid_driven_cavity
    echo 1 | ./NS_lid_driven_cavity                 # 1 = Upwind, 2 = Centred; writes EnsightOutput.case + *.vec/*.scl
    pvbatch make_video.py                           # -> frames/frame_0000.png ... (one per snapshot)
    ffmpeg -framerate 10 -i frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p velocity.mp4

The number of snapshots is `nstep/isto` (see `src/main.cpp`); the defaults `nstep = 20000`, `isto = 200`
give 100 frames. `pvbatch make_video.py <case-file> <frames-dir>` overrides the input/output paths.

## Contributing

Please read [CONTRIBUTING.md](https://github.com/DavidCico/Development-of-Navier-Stokes-code-by-projection-method/blob/master/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **David Cicoria** - *Initial work* - [DavidCico](https://github.com/DavidCico)

See also the list of [contributors](https://github.com/DavidCico/Development-of-Navier-Stokes-code-by-projection-method/graphs/contributors) who participated in this project.
