#
# Makefile
#
# snuppat
#

# Compiler.
FC = mpifort

# Flags.
 
# Optimization flags.
FFLAGS = -march=native -O3
MODFLAGS = -c -march=native -O3

# Profile flags.
ifeq ($(mode), p)
    FFLAGS = -march=native -O3 -g
    MODFLAGS = -c -march=native -O3 -g
endif

# Debug flags.
ifeq ($(mode), d)
    FFLAGS = -Og -ggdb3 -Wall -Wextra -Wconversion -pedantic -fbacktrace\
			 -ffpe-trap=invalid,zero,overflow -fbounds-check
    MODFLAGS = -c -Og -ggdb3 -Wall -Wextra -Wconversion -pedantic -fbacktrace\
			   -ffpe-trap=invalid,zero,overflow -fbounds-check
endif

# This Makefile
MAKEFILE = Makefile

# Executables.
EXE = snuppat

# Intermediate objects.
OBJECTS = loader integrator writer_reader mixer math_routines shell_operations

# Mods.
MODS = readvars_mod structures_mod integvars_mod

# All:
ALL = $(EXE) $(MODS:=.o) $(OBJECTS:=.o) tags

#--------------------------------------------

# All rule
all: $(ALL) $(MAKEFILE)

# Dependencies.
loader.o: readvars_mod.o structures_mod.o math_routines.o
integrator.o: integvars_mod.o structures_mod.o math_routines.o mixer.o
writer_reader.o: readvars_mod.o structures_mod.o shell_operations.o
shell_operations.o: readvars_mod.o structures_mod.o math_routines.o
mixer.o: structures_mod.o
math_routines.o: structures_mod.o
$(EXE): $(OBJECTS:=.o) $(MODS:=.o)

# Rules

# Ctags
tags: *.f90 *.py $(MAKEFILE)
	ctags --fortran-kinds=+L *.f90 *.py

# Executable compilation.
%: %.f90
	$(FC) $(FFLAGS) $(filter-out $(MAKEFILE), $^) -o $@

# Objects compilations.
%.o: %.f90
	$(FC) $(MODFLAGS) $<

# Cleaning.
clean:
	rm -f $(filter-out tags, $(ALL)) $(MODS:=.mod) $(OBJECTS:=.mod)
