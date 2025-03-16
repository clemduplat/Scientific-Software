----Question 1----

File :  -'program exponential_growth2.f95' : you must change precision= sp or dp based on the precision you want
Run ? : - gfortran -c program exponential_growth2.f95
Now you must fill in some variables and the terminal will write an answer based on your precision

----Question 2------

Here are the instructions for compiling my SIQRD programma:
The code contains:
	- Euler_F.f95 : contain the implementation of forward_euler, 
	- heun.f95  : contain the implementation of heun
	- back_euler.f95 : contain the implementation of backward_euler (forward and backward error)
	- Derivatives.f95: module containing the Deriviatives functions, Jacobian Function and identity matrix maker function
	- param.f95 : module containing the parameters of the program
To run it, you must have access to solver_gfortran.mod, solver_gfortran.o (or if you want to use another compiler, you must change the use solver_gfortran line 2 in SIQRD.f95 by the solver you need and proceed the same way)
Run ?:
	- gfortran -c param.f95
	- gfortran -c Derivatives.f95
	- gfortran -c "method_name".f95
	- gfortran -o method "method_name".o param.o Derivatives.o (solver_gfortran.o only for back_euler)
	- ./method
"method_name" is either Euler_F.f95, heun.f95 or back_euler.f95
Now you have access to a simulation_data.data file containing the evolution t S I Q R D at each grid point for the method
