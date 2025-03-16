#ifndef SOLVERS_HPP
#define SOLVERS_HPP

#include <iostream>
#include <fstream>
//#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <algorithm>
#include <cmath>
#include <functional>

using namespace boost::numeric::ublas;
//---Define Functors to make it general----//
using DerivativeFunction = std::function<vector<double>(const vector<double>&)>;
using JacobianFunction = std::function<matrix<double>(const vector<double>&)>;

//---Forward Euler solver----//
vector<double> solver_forward(const vector<double>& x, double dt,DerivativeFunction derivative) {
    vector<double> x_1(x.size()), fk(x.size());
    fk = derivative(x);
    x_1.assign(x + dt*fk);
    return x_1;
}


//-----Solver for Heun----/
vector<double> solver_heun(const vector<double>& x,double dt,DerivativeFunction derivative) {
    vector<double> x_1(x.size()), x_k(x.size()), fk(x.size()), fk_k(x.size());
    fk = derivative(x);
    x_k.assign(x+dt*fk);
    fk_k = derivative(x_k);
    x_1.assign(x + dt * 0.5 * (fk + fk_k));
    return x_1; 
}
//----Solver Backward-----//
vector<double> solver_backward(double dt, const vector<double>& x_0, DerivativeFunction derivative, JacobianFunction jacobian) {
    const double custom_tolerance = 1.0E-6;
    const int max_iterations = 100;
    vector<double> x = x_0;
    vector<double> x_new(x.size());
    vector<double> b(x.size());
    vector<double> f(x.size());
    matrix<double> J(x.size(), x.size());
    permutation_matrix<std::size_t> pm(x.size());
    for (int k = 0; k < max_iterations; ++k) {
        f = derivative(x);
        J = jacobian(x) * dt - identity_matrix<double>(x.size());
        b.assign(x_0 + dt * f - x);
        if (norm_2(b) < custom_tolerance * norm_2(x)) {
            break; 
        }
        lu_factorize(J, pm);
        lu_substitute(J, pm, b);
        x_new.assign(x - b); 
        if (norm_2(x_new - x) < custom_tolerance) {
            break; 
        }
        x.swap(x_new);
    }
    return x;
}



#endif