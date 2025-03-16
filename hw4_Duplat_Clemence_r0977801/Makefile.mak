# Variables
CC := g++ # Compiler
CFLAGS := -Wall -std=c++17 -O3 # Compiler flags


# Phony target for clean-up
.PHONY: clean
#SIQRDparams.hpp 
all: simulation1 simulation2 estimation1
# Compile simulation1.cpp to an object file
simulation1: simulation1.cpp SIQRDmodel.hpp solvers.hpp 
	$(CC) $(CFLAGS) -o simulation1 simulation1.cpp
simulation2: simulation2.cpp ODEmodel.hpp solvers.hpp 
	$(CC) $(CFLAGS) -o simulation2  simulation2.cpp

estimation: estimation1.cpp SIQRDmodel.hpp solvers.hpp gradient_approximation.hpp BFGS_algo.hpp
	$(CC) $(CFLAGS) -o estimation1  estimation1.cpp



# Clean-up command
clean:
	rm -f simulation1 simulation2 estimation1 *.o
