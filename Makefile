GF := gfortran
GFFLAGS := -o3 -std=f95

EXE_DIR    := bin/
OBJ_DIR    := obj/
SRC_DIR    := src/
EXECUTABLE := $(EXE_DIR)1D_solver
SRC_FILES  := $(wildcard src/M_precision.f90) \
							$(wildcard src/M_pin.f90) \
							$(wildcard src/M_domain.f90) \
							$(wildcard src/M_spatial_recon.f90) \
							$(wildcard src/M_integration.f90) \
							$(wildcard src/M_output.f90) \
							$(wildcard src/main.f90)

# Generally useful targets

.PHONY : all dirs clean

all : dirs $(EXECUTABLE)

dirs : $(EXE_DIR)

$(EXE_DIR):
	mkdir -p $(EXE_DIR)

$(EXECUTABLE) : 
	$(GF) $(GFFLAGS) -o $@ $(SRC_FILES)

clean :
	rm -rf $(EXECUTABLE)
	rm -rf *mod