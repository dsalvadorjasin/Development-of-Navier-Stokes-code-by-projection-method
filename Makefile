# Build of the lid-driven cavity solver (C++, sources in src/).
#
#   make                        -> ./NS_lid_driven_cavity
#   make check                  -> run both schemes and compare every output
#                                  with the Fortran reference oracle (reference/)
#   make clean
#
# The original Fortran program is kept in fortran/ as the reference oracle
# (see PORTING.md).  Two extra targets link C++ units into it to validate a
# port unit by unit:
#   make hybrid PORTED="timestep boundary_conditions"
#       -> ./NS_lid_driven_cavity_hybrid : Fortran program in which ONLY the
#          listed units are replaced by their C++ port (src/<unit>.cpp)
#   make check-hybrid PORTED="timestep"
#       -> builds the hybrid and compares both schemes with reference/
# Unit names are the source file stems (reading_data ... ensight_cavity).

CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
FC       = gfortran
FFLAGS   = -fdefault-real-8 -O3

SRC_DIR   = src
F90_DIR   = fortran
BUILD_DIR = build
EXE       = NS_lid_driven_cavity
HYBRID    = NS_lid_driven_cavity_hybrid

UNITS = reading_data mesh_creation initial_conditions boundary_conditions timestep \
        discretisation_upwind discretisation_centre calcul_2nd_membre matgen_cavity solver \
        calcul_pressure_field calcul_new_velocity velocities_output ensight_cavity
COMMON = variables support dump_fields

CPP_OBJS    = $(addprefix $(BUILD_DIR)/cpp/,$(addsuffix .o,main $(COMMON) $(UNITS)))
HYB_CPP_OBJS = $(addprefix $(BUILD_DIR)/hybrid/,$(addsuffix .o,$(COMMON) $(PORTED)))
FORTRAN_UNITS = $(filter-out $(PORTED),$(UNITS))
HYB_F_OBJS  = $(addprefix $(BUILD_DIR)/hybrid/f_,$(addsuffix .o,variables_module main $(FORTRAN_UNITS)))

.PHONY: all cpp hybrid check check-hybrid clean cleanall

all: cpp
cpp: $(EXE)

$(EXE): $(CPP_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD_DIR)/cpp/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/*.hpp | $(BUILD_DIR)/cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---- hybrid: Fortran program + selected C++ units ------------------------
hybrid: $(HYBRID)

$(HYBRID): $(HYB_F_OBJS) $(HYB_CPP_OBJS)
	$(CXX) $(CXXFLAGS) $^ -lgfortran -o $@

$(BUILD_DIR)/hybrid/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/*.hpp | $(BUILD_DIR)/hybrid
	$(CXX) $(CXXFLAGS) -DNS_HYBRID -c $< -o $@

# the Fortran module must be compiled first (it produces variables.mod)
$(BUILD_DIR)/hybrid/f_%.o: $(F90_DIR)/%.f90 $(BUILD_DIR)/hybrid/f_variables_module.o | $(BUILD_DIR)/hybrid
	$(FC) $(FFLAGS) -J$(BUILD_DIR)/hybrid -c $< -o $@

$(BUILD_DIR)/hybrid/f_variables_module.o: $(F90_DIR)/variables_module.f90 | $(BUILD_DIR)/hybrid
	$(FC) $(FFLAGS) -J$(BUILD_DIR)/hybrid -c $< -o $@

$(BUILD_DIR)/cpp $(BUILD_DIR)/hybrid:
	mkdir -p $@

# ---- validation against the reference oracle ------------------------------
check: $(EXE)
	tools/check_against_reference.sh ./$(EXE)

check-hybrid: $(HYBRID)
	tools/check_against_reference.sh ./$(HYBRID)

clean:
	rm -rf $(BUILD_DIR) $(EXE) $(HYBRID)

# also remove the simulation outputs written in the current directory
cleanall: clean
	rm -f *.case *.geo *.vec *.scl coordinates.txt u_velocity.txt v_velocity.txt omega_boundary.txt fields_*.txt
