#include "solvers.hpp"
#include "SIQRDmodel.hpp"
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <functional>
using namespace boost::numeric::ublas;
using DerivativeFunction = std::function<vector<double>(const vector<double>&)>;
using JacobianFunction = std::function<matrix<double>(const vector<double>&)>;


//----Compute the derivatives----//
int main(int argc, char const *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " N T\n";
        return 1;
    }
    SIQRDParams params;
    //r-----Extract parameters as asked-----//
    double N = std::stoi(argv[1]);
    double T = std::stod(argv[2]);
    double dt = T / static_cast<double>(N);
    std::ifstream paramFile("parameters.in");
    if (!paramFile.is_open()) {
        std::cerr << "Failed to open parameters.in\n";
        return 1;
    }
    double S0;
    double I0;
    paramFile >> params.beta >> params.mu >> params.gamma
              >> params.alpha >> params.delta
              >> S0 >> I0;
    paramFile.close();
    vector<double> stateForward(5), stateHeun(5), stateBackward(5);
    stateForward(0) = stateHeun(0) = stateBackward(0) = S0; 
    stateForward(1) = stateHeun(1) = stateBackward(1) = I0; 
    stateForward(2) = stateHeun(2) = stateBackward(2) = 0.0; 
    stateForward(3) = stateHeun(3) = stateBackward(3) = 0.0; 
    stateForward(4) = stateHeun(4) = stateBackward(4) = 0.0; 
    params.delta=0.0; 
    SIQRDModel model_f(params,stateForward);
    params.delta=0.2;
    SIQRDModel model_b(params,stateBackward);
    params.delta=0.9;
    SIQRDModel model_h(params,stateHeun);
    DerivativeFunction derivativeFunc_f = [&model_f](const vector<double>& state) {
        return model_f.calculate_derivatives(state);
    };
    DerivativeFunction derivativeFunc_h = [&model_h](const vector<double>& state) {
        return model_h.calculate_derivatives(state);
    };
    DerivativeFunction derivativeFunc_b = [&model_h](const vector<double>& state) {
        return model_h.calculate_derivatives(state);
    };
    //only need jacobian for backwardeuler
    JacobianFunction jacobian = [&model_b](const vector<double>& state) {
        return model_b.calculate_jacobian(state);
    };
    

    std::ofstream outputFileForward("fwe_no_measures.out");
    std::ofstream outputFileHeun("heun_lockdown.out");
    std::ofstream outputFileBackward("bwe_quarantine.out");
    
    if (!outputFileForward.is_open() || !outputFileHeun.is_open() || !outputFileBackward.is_open()) {
        std::cerr << "Failed to open files for writing\n";
        return 1;
    }

    for (int k = 0; k < N; ++k) {
        double t = k * dt;
    
        stateForward = solver_forward(stateForward, dt, derivativeFunc_f);

        stateHeun = solver_heun(stateHeun, dt, derivativeFunc_h);

        stateBackward = solver_backward(dt,stateBackward, derivativeFunc_b, jacobian);

        outputFileForward << t << " " << stateForward(0)<< " " << stateForward(1)<< " " << stateForward(2) << " " << stateForward(3)<< " " << stateForward(4) << std::endl;
        outputFileHeun << t << " " << stateHeun(0) << " " << stateHeun(1) << " " << stateHeun(2) << " " << stateHeun(3) << " " << stateHeun(4) <<std::endl;
        outputFileBackward << t << " " << stateBackward(0)<< " " << stateBackward(1)<< " " << stateBackward(2) << " " << stateBackward(3) << " " << stateBackward(4) <<std::endl;
    }

    
    double sumForward = sum(stateForward);
    double sumHeun = sum(stateHeun);
    double sumBackward = sum(stateBackward);

    std::cout << "Sum of the last elements (Forward Euler): " << sumForward << std::endl;
    std::cout << "Sum of the last elements (Heun's Method): " << sumHeun << std::endl;
    std::cout << "Sum of the last elements (Backward Euler): " << sumBackward << std::endl;

    outputFileForward.close();
    outputFileHeun.close();
    outputFileBackward.close();
    
    return 0;
}
