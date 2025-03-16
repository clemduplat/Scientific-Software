#-------------------------------------Clémence Duplat r0977801--------------------------#
#------------------README-------------#
Welcome in the SIQRD simulation model
This zip contains the following modules
	- param.f95 : initialize the parameters need
	- Derivatives.f95 : calculate the derivatives/Jacobian of the SIQRD model
	- solver.f95 : generic interface for solver function that uses LAPACK library (-llapack)
	- solvers.f95 : contains the implementation of Forward_Euler, Heun and Backward_Euler method

The zip contains the following main programs 
	- Question2.f95 : uses of the method with the parameters given in the question and N and T asked to the user
				the choice of method is also asked to the user
	- tweaking.f90: optimisation problem that find an optimal beta based on a range of initial Delta(beta) and an initial target
				the choice of method is also asked to the user
	-satbility.f90: find the real part of the eigenvalues of the Jacobian to analyse stability 
The zip contains a Makefile:
	- You can change the program you want to launch (already put in comment)
	- You can change the flags you want to use
	- You can change the compiler to gfortran, ifort or nagfor
	- make: to launch it
	- make clean: to delete the .o extensions