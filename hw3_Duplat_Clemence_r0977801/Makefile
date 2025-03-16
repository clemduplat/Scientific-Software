# Define the compiler
FC = ifort #ifort #nagfor

# Define flags for the compiler
#FFLAGS = -O3
#FFLAGS = -O0
FFLAGS = -Og -fbounds-check
#FFLAGS = -O3 -funroll-loops
#FFLAGS = -Ofast -march=native -ffast-math
# Define the target executable
TARGET = session_ifort

#TARGET = i2_block
#TARGET = blas

# Define the main source file
MAIN_SRC = mm_driver.f90
#MAIN_SRC = i2.f90
#MAIN_SRC=blas_divide_and_conquer.f90
# Define additional source files
ADDITIONAL_SRC = matrixop.f90 timings.f90

# Define object files (compiled additional source files)
OBJ = $(ADDITIONAL_SRC:.f90=.o)

# Rule for creating the executable
$(TARGET): $(MAIN_SRC) $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(MAIN_SRC) $(OBJ) -lblas

# Rule for creating object files
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Clean rule
clean:
	rm -f $(TARGET) $(OBJ)

# Phony targets
.PHONY: clean
