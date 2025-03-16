# Define the compiler
FC = gfortran
#FC = ifort

# Define compiler flags
FFLAGS = -O0 
#FFLAGS = -Og
#FFLAGS = -O2
#FFLAGS = -O3
#FFLAGS = -Ofast
# Define the LAPACK library flags
LAPACK_LIBS = -llapack

# Object files
OBJS = tweaking.o timing.o solvers.o solver.o Derivatives.o param.o
#OBJS = stability.o solvers.o solver.o Derivatives.o param.o
# Executable name
EXEC = tweaking
#EXEC = stability

# Default target
all: $(EXEC)

# Linking the executable
$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LAPACK_LIBS)

# Compiling the source files
param.o: param.f95
	$(FC) $(FFLAGS) -c $<

Derivatives.o: Derivatives.f95
	$(FC) $(FFLAGS) -c $<

solver.o: solver.f95
	$(FC) $(FFLAGS) -c $<

solvers.o: solvers.f95
	$(FC) $(FFLAGS) -c $<

timing.o: timing.f90
	$(FC) $(FFLAGS) -c $<
#stability.o: stability.f90
tweaking.o: tweaking.f90
	$(FC) $(FFLAGS) -c $<

# Clean target
clean:
	rm -f $(OBJS) $(EXEC)

# Mark all object files as phony
.PHONY: all clean $(OBJS)
