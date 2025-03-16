#----README----#
This is the assignement 4 in C++ for the Scientific Software Course at KULeuven
In this zip you can find:
    - SIQRDmodel.hpp: class definition of SIQRD model
    - EDOmodel.hpp: class with the general EDO given in the assignement
    - solvers.hpp: IVP for general EDO that takes a specific derivative and Jacobian based on the model
    - simulation1.cpp: simulation of the IVPs with SIQRDmodel. Note that we still have to change the delta manually, I hadn't the time to finish my Object Oriented SIQRD model
    - simulation2.cpp: simulation of the IVPs with EDOmodel.
    - gradient_approximation.hpp: contains the LSE function calculation, gradient of the LSE function calculation, linesearch function that returns optimal step size            and finally the BFGS algorithm based on all those functions
    - estimation1.cppp: applied BFGS algorithm to an observation file. Doesn't work as intended
    - Makefile: make all, make simulation1, make simulation2, make estimation1
    Note that:  ./simulation N T
    - .in files used for the simulations and estimations
    - report about the code
