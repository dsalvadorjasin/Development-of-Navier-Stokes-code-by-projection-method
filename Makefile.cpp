# C++ build of the lid-driven cavity solver (port of the Fortran code).
#
#   make -f Makefile.cpp                 -> ./NS_lid_driven_cavity_cpp   (pure C++)
#   make -f Makefile.cpp hybrid PORTED="timestep boundary_conditions"
#                                        -> ./NS_lid_driven_cavity_hybrid
#      Fortran program in which ONLY the listed units are replaced by their
#      C++ port (src/<unit>.cpp); all other units stay Fortran.  Used to
#      validate each ported unit in isolation against the reference oracle:
#          echo 1 | NS_DUMP_FIELDS=1 ./NS_lid_driven_cavity_hybrid
#          tools/compare_reference.py . reference/upwind
#   make -f Makefile.cpp check-hybrid PORTED="timestep"
#      builds the hybrid and runs both schemes against reference/{upwind,centred}
#   make -f Makefile.cpp check          -> same for the pure C++ executable
#
# Unit names are the Fortran file stems: reading_data mesh_creation
# initial_conditions boundary_conditions timestep discretisation_upwind
# discretisation_centre calcul_2nd_membre matgen_cavity solver
# calcul_pressure_field calcul_new_velocity velocities_output ensight_cavity

CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
FC       = gfortran
FFLAGS   = -fdefault-real-8 -O3

SRC_DIR   = src
BUILD_DIR = build
EXE       = NS_lid_driven_cavity_cpp
HYBRID    = NS_lid_driven_cavity_hybrid

UNITS = reading_data mesh_creation initial_conditions boundary_conditions timestep \
        discretisation_upwind discretisation_centre calcul_2nd_membre matgen_cavity solver \
        calcul_pressure_field calcul_new_velocity velocities_output ensight_cavity
COMMON = variables support dump_fields

CPP_OBJS    = $(addprefix $(BUILD_DIR)/cpp/,$(addsuffix .o,main $(COMMON) $(UNITS)))
HYB_CPP_OBJS = $(addprefix $(BUILD_DIR)/hybrid/,$(addsuffix .o,$(COMMON) $(PORTED)))
FORTRAN_UNITS = $(filter-out $(PORTED),$(UNITS))
HYB_F_OBJS  = $(addprefix $(BUILD_DIR)/hybrid/f_,$(addsuffix .o,variables_module main $(FORTRAN_UNITS)))

.PHONY: all cpp hybrid check check-hybrid clean

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
$(BUILD_DIR)/hybrid/f_%.o: %.f90 $(BUILD_DIR)/hybrid/f_variables_module.o | $(BUILD_DIR)/hybrid
	$(FC) $(FFLAGS) -J$(BUILD_DIR)/hybrid -c $< -o $@

$(BUILD_DIR)/hybrid/f_variables_module.o: variables_module.f90 | $(BUILD_DIR)/hybrid
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
