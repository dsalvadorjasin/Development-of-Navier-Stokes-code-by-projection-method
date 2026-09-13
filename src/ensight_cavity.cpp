// C++ port of ensight_cavity.f90: Ensight Gold binary output
// (write_result_ensight, WriteRectEnsightGeo, WriteEnsightVar, EnsightCase).
//
// The .geo/.vec/.scl files are Fortran `unformatted, access='direct'` files
// holding ONE record of length `reclength`: no record markers, the payload in
// write order (80-byte space-padded strings, 4-byte ints, 4-byte floats),
// padded with NUL bytes up to `reclength`.  Written here with a byte buffer
// that mirrors the items of the Fortran `write(unit,rec=1) ...` statement.
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "routines.hpp"
#include "variables.hpp"

namespace ns {

namespace {

// character(len=80): value left-justified, padded with blanks to 80 bytes
std::string fstring80(const std::string& s) {
    std::string out(80, ' ');
    out.replace(0, s.size() < 80 ? s.size() : 80, s, 0, 80);
    return out;
}

// one Fortran direct-access record: items are appended in write order
class DirectRecord {
public:
    void put_string80(const std::string& s) {
        const std::string padded = fstring80(s);
        buf_.insert(buf_.end(), padded.begin(), padded.end());
    }
    void put_int4(int v) {  // default integer, little-endian on x86-64
        const std::int32_t i = v;
        const unsigned char* b = reinterpret_cast<const unsigned char*>(&i);
        buf_.insert(buf_.end(), b, b + 4);
    }
    void put_real4(double v) {  // real(sngl(x),4)
        const float f = static_cast<float>(v);
        const unsigned char* b = reinterpret_cast<const unsigned char*>(&f);
        buf_.insert(buf_.end(), b, b + 4);
    }
    // open(..., access='direct', recl=reclength) + write(rec=1): the record is
    // padded with NUL bytes to recl (gfortran behaviour); the file is created
    // fresh / truncated.
    void write_file(const std::string& filename, int reclength) const {
        std::vector<unsigned char> rec(buf_);
        if (static_cast<int>(rec.size()) > reclength) {
            std::fprintf(stderr, "error: Ensight record of %zu bytes exceeds recl=%d (%s)\n",
                         rec.size(), reclength, filename.c_str());
            std::abort();
        }
        rec.resize(reclength, 0);
        std::ofstream f(filename, std::ios::binary | std::ios::trunc);
        f.write(reinterpret_cast<const char*>(rec.data()), rec.size());
    }

private:
    std::vector<unsigned char> buf_;
};

// Fortran Iw edit descriptor: right-justified in w columns, all '*' on overflow
std::string fint(int w, int v) {
    char tmp[32];
    std::snprintf(tmp, sizeof tmp, "%d", v);
    std::string s(tmp);
    if (static_cast<int>(s.size()) > w) return std::string(w, '*');
    return std::string(w - s.size(), ' ') + s;
}

// Fortran trim(): remove trailing blanks
std::string ftrim(const std::string& s) {
    const std::string::size_type e = s.find_last_not_of(' ');
    return e == std::string::npos ? std::string() : s.substr(0, e + 1);
}

}  // namespace

//    ******************************************************************************
//                      subrout  WriteRectEnsightGeo
//    ******************************************************************************
//    writes mesh data in Ensight's binary format for rectilinear geometry
void WriteRectEnsightGeo(int imin, int imax, int jmin, int jmax, int kmin, int kmax,
                         FArray1<double>& x1, FArray1<double>& x2, FArray1<double>& x3,
                         const std::string& FileName) {
    std::string binary_form, file_description1, file_description2, node_id, element_id, part,
        description_part, blocck;
    int i, j, k, npart, isize, jsize, ksize;
    int reclength;

    binary_form = "C Binary";
    file_description1 = "Ensight Model Geometry File Created by ";
    file_description2 = "WriteRectEnsightGeo Routine";
    node_id = "node id off";
    element_id = "element id off";
    part = "part";
    npart = 1;
    description_part = "Sortie Ensight";
    blocck = "block rectilinear";
    isize = (imax - imin + 1);
    jsize = (jmax - jmin + 1);
    ksize = (kmax - kmin + 1);

    reclength = 80 * 8 + 4 * (4 + isize + jsize + ksize);

    DirectRecord rec;
    rec.put_string80(binary_form);
    rec.put_string80(file_description1);
    rec.put_string80(file_description2);
    rec.put_string80(node_id);
    rec.put_string80(element_id);
    rec.put_string80(part);
    rec.put_int4(npart);
    rec.put_string80(description_part);
    rec.put_string80(blocck);
    rec.put_int4(isize);
    rec.put_int4(jsize);
    rec.put_int4(ksize);
    for (i = imin; i <= imax; ++i) rec.put_real4(x1(i));
    for (j = jmin; j <= jmax; ++j) rec.put_real4(x2(j));
    for (k = kmin; k <= kmax; ++k) rec.put_real4(x3(k));

    rec.write_file(ftrim(FileName) + ".geo", reclength);
}

//    ******************************************************************************
//    WriteEnsightVar writes result data in Ensight's format
//
//    ndv.........: number of dimension of the variable (1=>scalar   3=>vector)
//    var.........: data to be written, Fortran array var(ndv,nx,ny,nz)
//                  (column-major: element (m,i,j,k) at
//                   (m-1) + ndv*((i-1) + nx*((j-1) + ny*(k-1))))
//    Varname.....: word used to build filenames
//    imin,imax...: range of writting data in the x1 direction
//    jmin,jmax...: range of writting data in the x2 direction
//    kmin,kmax...: range of writting data in the x3 direction
//    ******************************************************************************
void WriteEnsightVar(int ndv, int nx, int ny, int nz, const std::vector<double>& var,
                     const std::string& Varname, int imin, int imax, int jmin, int jmax, int kmin,
                     int kmax) {
    std::string VarFileName;
    std::string part, blocck;
    int i, j, k, npart, m;
    int reclength;
    (void)nz;

    part = "part";
    npart = 1;
    blocck = "block rectilinear";

    if (ndv == 1) VarFileName = ftrim(Varname) + ".scl";
    if (ndv == 3) VarFileName = ftrim(Varname) + ".vec";

    reclength = 80 * 3 + 4 * (1 + (imax - imin + 1) * (jmax - jmin + 1) * (kmax - kmin + 1) * ndv);

    DirectRecord rec;
    rec.put_string80(VarFileName);
    rec.put_string80(part);
    rec.put_int4(npart);
    rec.put_string80(blocck);
    // ((((var(m,i,j,k), i=imin,imax), j=jmin,jmax), k=kmin,kmax), m=1,ndv)
    for (m = 1; m <= ndv; ++m)
        for (k = kmin; k <= kmax; ++k)
            for (j = jmin; j <= jmax; ++j)
                for (i = imin; i <= imax; ++i)
                    rec.put_real4(var[(m - 1) + ndv * ((i - 1) + nx * ((j - 1) + ny * (k - 1)))]);

    rec.write_file(VarFileName, reclength);
}

//    ******************************************************************************
//   EnsightCase helps to write a Ensight's case file
//
//    VarName.....: Name of the variable
//    GeoName.....: Name of the geometrie
//    VarType.....: 1 => Scalar       2 => Vector
//    nfile.......: number of result files (time steps)
//    ******************************************************************************
void EnsightCase(int nbVar, const std::vector<std::string>& VarName, const std::string& GeoName,
                 const std::vector<int>& VarType, int istep, int nstep, int isto) {
    int i, j, nfile;
    (void)istep;

    nfile = nstep / isto;

    std::ofstream f(ftrim(GeoName) + ".case", std::ios::trunc);

    // 10 format('FORMAT',/,'type: ensight gold',//,'GEOMETRY',/,'model:<TAB>',A,//,'VARIABLE')
    f << "FORMAT\n"
      << "type: ensight gold\n"
      << "\n"
      << "GEOMETRY\n"
      << "model:\t" << ftrim(GeoName) + ".geo" << "\n"
      << "\n"
      << "VARIABLE\n";

    // 15 format('scalar per node: ',A,'   ', A)
    // 25 format('vector per node: ',A,'   ', A)
    for (i = 1; i <= nbVar; ++i) {
        const std::string name = ftrim(VarName[i - 1]);
        if (nfile == 1) {
            if (VarType[i - 1] == 1) f << "scalar per node: " << name << "   " << name + ".scl" << "\n";
            if (VarType[i - 1] == 2) f << "vector per node: " << name << "   " << name + ".vec" << "\n";
        } else {
            if (VarType[i - 1] == 1)
                f << "scalar per node: " << name << "   " << name + "******.scl" << "\n";
            if (VarType[i - 1] == 2)
                f << "vector per node: " << name << "   " << name + "******.vec" << "\n";
        }
    }

    if (nfile > 1) {
        // 45 format(/,'TIME            ',/,'time set: 1     ',/,'number of steps:',i4,/,
        //           'filename start number:',i10/,'filename increment:',i10/,'time values: ')
        f << "\n"
          << "TIME            \n"
          << "time set: 1     \n"
          << "number of steps:" << fint(4, nfile) << "\n"
          << "filename start number:" << fint(10, isto) << "\n"
          << "filename increment:" << fint(10, isto) << "\n"
          << "time values: \n";
        // write(FileUnit,'(i5)') (j*isto,j=1,nfile)   -> one value per record
        for (j = 1; j <= nfile; ++j) f << fint(5, j * isto) << "\n";
    }
}

//-------------------------------------------------!
// Subroutine to write results in Ensight format   !
//-------------------------------------------------!
void write_result_ensight(FArray1<double>& xEnsight, FArray1<double>& yEnsight,
                          FArray2<double>& u, FArray2<double>& v, FArray2<double>& omega,
                          FArray2<double>& press, int nx, int ny, int nz, int istep, int isto,
                          int nstep) {
    constexpr int nbVar = 3;
    std::vector<double> velocity(static_cast<std::size_t>(3) * nx * ny * nz);  // (3,nx,ny,nz)
    std::vector<double> scalar(static_cast<std::size_t>(1) * nx * ny * nz);    // (1,nx,ny,nz)
    int nbImages;
    std::string caseName, speedName, omegaName, presName;
    std::vector<std::string> varName(nbVar);
    std::vector<int> varType(nbVar);
    double zEnsight_store[1];
    FArray1<double> zEnsight(zEnsight_store, 1, 1);
    int i, j, imin, imax, jmin, jmax, kmin, kmax;
    char step6[16];

    // velocity(m,i,j,1) / scalar(1,i,j,1) in the column-major layout of the Fortran locals
    auto vel = [&](int m, int ii, int jj) -> double& {
        return velocity[(m - 1) + 3 * ((ii - 1) + nx * (jj - 1))];
    };
    auto sca = [&](int ii, int jj) -> double& { return scalar[(ii - 1) + nx * (jj - 1)]; };

    caseName = "EnsightOutput";
    varName[0] = "velocity";
    varType[0] = 2;
    varName[1] = "omega";
    varType[1] = 1;
    varName[2] = "press";
    varType[2] = 1;
    nbImages = nstep / isto;
    (void)nbImages;

    zEnsight.fill(0.);

    kmin = 1;
    kmax = 1;
    jmin = 1;
    jmax = ny;
    imin = 1;
    imax = nx;
    if (nstep / isto == 1) {
        speedName = "velocity";
        omegaName = "omega";
        presName = "press";
    } else {
        std::snprintf(step6, sizeof step6, "%06d", istep);  // '(A,I6.6)'
        speedName = std::string("velocity") + step6;
        omegaName = std::string("omega") + step6;
        presName = std::string("press") + step6;
    }
    EnsightCase(nbVar, varName, caseName, varType, istep, nstep, isto);
    WriteRectEnsightGeo(imin, imax, jmin, jmax, kmin, kmax, xEnsight, yEnsight, zEnsight, caseName);

    for (j = 1; j <= ny; ++j) {
        for (i = 1; i <= nx; ++i) {
            vel(1, i, j) = u(i, j);
            vel(2, i, j) = v(i, j);
            vel(3, i, j) = 0.;
        }
    }

    WriteEnsightVar(3, nx, ny, nz, velocity, speedName, imin, imax, jmin, jmax, kmin, kmax);
    for (j = 1; j <= ny; ++j) {
        for (i = 1; i <= nx; ++i) {
            sca(i, j) = omega(i, j);
        }
    }

    // the Fortran passes jmin where kmin is expected (both are 1); kept as is
    WriteEnsightVar(1, nx, ny, nz, scalar, omegaName, imin, imax, jmin, jmax, jmin, kmax);
    for (j = 1; j <= ny; ++j) {
        for (i = 1; i <= nx; ++i) {
            sca(i, j) = press(i, j);
        }
    }
    WriteEnsightVar(1, nx, ny, nz, scalar, presName, imin, imax, jmin, jmax, kmin, kmax);
}

}  // namespace ns

#ifdef NS_HYBRID
// entry point used when this routine replaces ensight_cavity.f90 inside the Fortran program
extern "C" void write_result_ensight_(double* x, double* y, double* u, double* v, double* omega,
                                      double* press, int* nx, int* ny, int* nz, int* istep,
                                      int* isto, int* nstep) {
    ns::FArray1<double> x_v(x, 1, *nx), y_v(y, 1, *ny);
    ns::FArray2<double> u_v(u, 1, *nx, 1, *ny), v_v(v, 1, *nx, 1, *ny),
        omega_v(omega, 1, *nx, 1, *ny), press_v(press, 1, *nx, 1, *ny);
    ns::write_result_ensight(x_v, y_v, u_v, v_v, omega_v, press_v, *nx, *ny, *nz, *istep, *isto,
                             *nstep);
}
#endif
