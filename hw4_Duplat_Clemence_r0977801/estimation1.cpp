

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include "solvers.hpp"
#include "SIQRDmodel.hpp"
#include "gradient_approximation.hpp"


using namespace boost::numeric::ublas;

std::pair<int,matrix<double>> read_observations(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }
    

    int num_observations = 0, num_variables = 5; // Assuming 5 variables per observation
    file >> num_observations; 

    matrix<double> observations(num_observations, num_variables);
    for (int i = 0; i < num_observations; ++i) {
        for (int j = 0; j < num_variables; ++j) {
            if (!(file >> observations(i, j))) {
                throw std::runtime_error("Error reading data at row " + std::to_string(i) + " column " + std::to_string(j));
            }
        }
    }
    
    return make_pair(num_observations,observations);
}





int main(int argc, char const *argv[]) {
    
    std::string filename = "observations1.in"; 
    //std::string filename = "observations2.in"; 
    std::pair<int, matrix<double>> result = read_observations(filename);

    int num_observations = result.first;
    matrix<double> observations = result.second;
    vector<double> initial_state = row(observations,0);
    vector<double> initial_params(5);
    initial_params(0)=0.32;//0.32
    initial_params(1)=0.03;
    initial_params(2)= 0.151;
    initial_params(3) =0.004;
    initial_params(4)= 0.052;
    SIQRDParams siqrdParams = {initial_params[0], initial_params[1], initial_params[2], initial_params[3], initial_params[4]};
    
    SIQRDModel model(siqrdParams, initial_state);
    model.set_initial_state(initial_state);
    
    auto modelFunction = [&model, &siqrdParams, &num_observations](const vector<double>& params) -> matrix<double> {
        //model.set_params({params[0], params[1], params[2], params[3], params[4]}); 
        
        siqrdParams.beta=params[0];
        siqrdParams.mu=params[1];
        siqrdParams.gamma=params[2];
        siqrdParams.alpha=params[3];
        siqrdParams.delta=params[4];
        //std :: cout << siqrdParams.alpha; // change a little bit but not too much
        //SIQRDModel model(siqrdParams, model.get_initial_state());
        vector<double> state = model.get_initial_state();
        int timeSteps = num_observations; // Adjust as necessary
        const double dt = 1/8; 

        matrix<double> predictions(timeSteps, state.size());
        DerivativeFunction derivativeFunc_f = [&model](const vector<double>& state) {
            return model.calculate_derivatives(state);
        };
        for (int step = 0; step < timeSteps; ++step) {
            state = solver_heun(state, dt, derivativeFunc_f);
            
            row(predictions, step) = state;
        }
        //std :: cout << predictions;
        return predictions;
    };

    int T=201;
    vector<double> optimized_params = BFGS(observations, modelFunction,initial_params, T,100);
    std::cout << "Optimized Parameters: ";
    for (std::size_t i = 0; i < optimized_params.size(); ++i) {
        std::cout << optimized_params(i) << " ";
    }
    std::cout << std::endl;
    // Output results and save predictions.
    
    return 0;
}