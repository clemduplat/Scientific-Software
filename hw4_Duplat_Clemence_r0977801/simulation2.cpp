#include "solvers.hpp"
#include "ODEmodel.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <numeric>
#include <functional>
#include <iterator>
#include <algorithm>

using namespace boost::numeric::ublas;
using DerivativeFunction = std::function<vector<double>(const vector<double>&)>;
using JacobianFunction = std::function<matrix<double>(const vector<double>&)>;



int main(int argc, char const *argv[]) {
    // Initialize state vector with initial values using iota
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " N T\n";
        return 1;
    }
    //r-----Extract parameters as asked-----//
    double N = std::stoi(argv[1]);
    double T = std::stod(argv[2]);
    double dt = T / static_cast<double>(N);
    vector<double> x_f(50),x_h(50),x_b(50);
    std::iota(std::begin(x_f), std::end(x_f), 1.0); 
    std::transform(std::begin(x_f), std::end(x_f), std::begin(x_f), [](double val) { return val * 0.01; });
    x_h = x_f; 
    x_b = x_f; 
    ODEModel model(x_f.size());
    DerivativeFunction derivativeFunc = [&model](const vector<double>& state) -> vector<double> {
        return model.calculateSystemDerivatives(state);
    };

    // Lambda function for jacobian
    JacobianFunction jacobian = [&model](const vector<double>& state) -> matrix<double> {
        return model.calculateJacobian(state);
    };
    
    std::ofstream fwe_file("fwe_simulation2.out");
    std::ofstream heun_file("heun_simulation2.out");
    std::ofstream bwe_file("bwe_simulation2.out");
    
    for (int i = 0; i <= N; ++i) {

        if (i > 0) {
            x_f = solver_forward(x_f,dt,derivativeFunc);
            x_h = solver_heun(x_h,dt,derivativeFunc);
            x_b= solver_backward(dt, x_b, derivativeFunc,jacobian);
        }
       

        
        fwe_file << i*dt ;
        heun_file << i * dt;
        bwe_file << i*dt ;
        for (size_t j = 0; j < x_f.size(); ++j) {
            fwe_file << " " << x_f(j);
            heun_file << " " << x_h(j);
            bwe_file << " " << x_b(j);

        }
        fwe_file << std::endl;
        heun_file << std::endl;
        bwe_file << std::endl;
    }

    // Close output files
    //double sumForward = sum(x_f);
    //double sumHeun = sum(x_h);
    //double sumBackward = sum(x_b);

    //std::cout << "Sum of the last elements (Forward Euler): " << sumForward << std::endl;
    //std::cout << "Sum of the last elements (Heun's Method): " << sumHeun << std::endl;
    //std::cout << "Sum of the last elements (Backward Euler): " << sumBackward << std::endl;
    fwe_file.close();
    heun_file.close();
    bwe_file.close();

    return 0;
}
