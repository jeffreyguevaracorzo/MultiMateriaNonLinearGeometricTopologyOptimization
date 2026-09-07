FC      = gfortran
MOD_DIR = $(BUILD_DIR)/mod
FFLAGS  = -O3 -fbounds-check -fbacktrace -fcheck=all -g -Wall -Wextra -Wrealloc-lhs-all -fopenmp -J$(MOD_DIR) -I$(MOD_DIR)

# Library paths.
#  - macOS + Homebrew: uncomment the LDFLAGS block below and adjust the versions.
#  - Linux with system packages: leave LDFLAGS empty.

LDFLAGS = -L/opt/homebrew/Cellar/metis/5.1.0/lib \
        -L/opt/homebrew/Cellar/libomp/22.1.2/lib
LIBS    = -llapack -lblas -lmetis -fopenmp

SRC_DIR     = src
BUILD_DIR   = build
BIN_DIR     = bin
TARGET      = $(BIN_DIR)/programa

# Compilation order matters (Fortran module dependencies)
SRC_FILES = \
	src/solver/sdeps90.f90 \
	src/solver/common.f \
	src/solver/common90.f90 \
	src/solver/hsl_ma86s.f90 \
	src/solver/hsl_ma86d.f90 \
	src/mma/MMA_Routines.f90 \
	src/mma/MMA_Interface.f90 \
	src/solver/Solver_MA86Module.f90 \
	src/base/Base_Module.f90 \
	src/fea/Base_FEA_MMNL_Module.f90 \
	src/fea/FEA_MMNL_Module.f90 \
	src/optimization/MMNL_Optimization_Module.f90 \
	src/postprocessing/MMNL_Paraview_Module.f90 \
	src/main/MainMultiMaterialNL.f90

OBJ_FILES = $(patsubst $(SRC_DIR)/%.f90, $(BUILD_DIR)/%.o, $(SRC_FILES))

# NOTE: do NOT use "make -j". The bundled HSL sources define modules consumed by
# sibling files in the same directory, and parallel compilation races on the .mod
# files. Serial "make" is correct and fast enough.
.NOTPARALLEL:

all: $(TARGET)

$(TARGET): $(OBJ_FILES) | $(BIN_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR) $(MOD_DIR)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(MOD_DIR):
	mkdir -p $(MOD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)/* *.mod \
	       output/paraview/*.case output/paraview/*.geom output/paraview/*.esca

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run
