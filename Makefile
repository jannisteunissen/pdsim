F90 := gfortran
FFLAGS := -O2 -g -std=f2008 -fopenmp -Wall -Wextra -J src -fcheck=all
PROGS := pdsim
LIBS := interp_unstructured particle_core
LIBDIRS := interpolate_unstructured particle_core
INCDIRS := interpolate_unstructured particle_core

.PHONY:	all clean

all: 	$(PROGS)

clean:
	$(RM) $(progs) src/*.o src/*.mod

# Dependency information
$(PROGS): src/m_config.o src/m_photoi.o src/m_pdsim.o

# How to get .o object files from .f90 source files
src/%.o: src/%.f90
	$(F90) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
%: src/%.o
	$(F90) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
