F90 := gfortran
FFLAGS := -O2 -g -std=f2008 -fopenmp -Wall -Wextra -J build -fcheck=all
LIBS := interp_unstructured particle_core
LIBDIRS := interpolate_unstructured particle_core
INCDIRS := interpolate_unstructured particle_core

.PHONY:	all clean

all: 	pdsim

clean:
	$(RM) build/*.o build/*.mod build/pdsim pdsim

# How to get .o object files from .f90 source files
build/%.o: src/%.f90
	$(F90) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

# How to get executables from .o object files
build/%: build/%.o
	$(F90) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

# Copy step for executables
%: build/%
	cp $^ $@

# Dependency information generated with
# fortdepend -f src/*.f90 -i m_cross_sec m_gas m_interp_unstructured m_particle_core m_random m_units_constants omp_lib -b build -w -o dependencies.make
include dependencies.make
