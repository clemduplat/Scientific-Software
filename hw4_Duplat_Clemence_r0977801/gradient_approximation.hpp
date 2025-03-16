#ifndef GRADIENT_APPROXIMATION_HPP
#define GRADIENT_APPROXIMATION_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <string>

using namespace boost::numeric::ublas;
//---General LSE----//
double LSE(const matrix<double>& observations, const matrix<double>& predictions, int T) {
    double lse = 0.0;
    for (int i = 0; i < T; ++i) {
        vector<double> error = row(observations, i) - row(predictions, i);
        double population = std::accumulate(row(observations, i).begin(), row(observations, i).end(), 0.0);
        lse += std::pow(norm_2(error),2) / std::pow(population,2);
    }
    return lse / T;
}
//---Gradient of my LSE---//
vector<double> gradient_of_LSE(const matrix<double>& observed, 
                                               std::function<matrix<double>(const vector<double>&)> model, 
                                               const vector<double>& parameters, 
                                               double epsilon,int T) {
    vector<double> gradient(parameters.size());
    vector<double> params_plus = parameters;
    double error = LSE(observed, model(parameters),T);
    for (unsigned i = 0; i < parameters.size(); ++i) {
        params_plus(i) += epsilon;
        double error_plus = LSE(observed, model(params_plus),T);
        gradient(i) = (error_plus - error) / (params_plus(i)* epsilon);
    }
    return gradient;
}
//----Inexact lineSearch----//
double lineSearch(const matrix<double>& observed, 
                  std::function<matrix<double>(const vector<double>&)> model, 
                  vector<double> parameters, 
                  double alpha_init, double c1, double c2, int T,vector<double> d) {
    double alpha = alpha_init;
    vector<double> gradient = gradient_of_LSE(observed, model, parameters, 1e-5, T);
    vector<double> new_parameters;
    double lse_old = LSE(observed, model(parameters), T);
    for (int iter = 0; iter < 1000; ++iter) {
        new_parameters = parameters + alpha * d;
        double lse_new = LSE(observed, model(new_parameters), T);
        //std :: cout <<lse_new; // doesn't change weiiiiird
        // Armijo condition
        if (lse_new < lse_old + c1 * alpha * inner_prod(d, gradient)) {
            std :: cout << "Hello I met Armajio in line search"; // not went here probleù
    
            break;
        }

    
        alpha=alpha/2;
    }

    return alpha;
}


vector<double> BFGS(const matrix<double>& observed,
                    std::function<matrix<double>(const vector<double>&)> model,
                    vector<double> parameters,
                    int T, int max_iterations, double tolerance = 1e-10) {
    
    const size_t n = parameters.size();
    // Declaration here to avoid useless storage---//
    matrix<double> H = identity_matrix<double>(n);  
    vector<double> grad = gradient_of_LSE(observed, model, parameters, 1e-5, T);
    vector<double> new_parameters(n);
    vector<double> new_grad(n);
    vector<double> d(n), s(n), y(n);
    matrix<double> temp_matrix(n, n);
    matrix<double> s_y_outer;
    double y_s_inner;
    int iter;
    double alpha =1.0;
    double c1=1e-4;
    double c2=0.9;
    double epsiolon=1e-5;

    for (iter = 0; iter < max_iterations; ++iter) {
        noalias(d) = -prod(H, grad); 
        alpha = lineSearch(observed, model, parameters, alpha, c1, c2, T,d);
        //std :: cout << alpha;
        
        new_parameters.assign(parameters +  alpha * d);
        if (alpha*norm_2(d)/norm_2(new_parameters)) {
            std :: cout <<"Hello I want in the stopping criteria in iter= ";
            std :: cout << iter ;
            std :: cout <<"  ";
            // Direct in the stopping criteria not normal !!!!!
            break; 
        }
        s.assign(alpha * d);
        //noalias(s) = new_parameters - parameters;  
        new_grad = gradient_of_LSE(observed, model, new_parameters, epsiolon, T);
        noalias(y) = new_grad - grad;
        
        

        //----Updating Hessian matrix---//
        y_s_inner = inner_prod(y, s);
        s_y_outer = outer_prod(s, y);
        temp_matrix.assign(prod(H, s_y_outer));
        H -= prod(temp_matrix, H) / y_s_inner;
        H += outer_prod(y, y) / y_s_inner;
        //std :: cout << new_grad;
        noalias(grad) = new_grad;
        //std:: cout << grad;
        parameters.assign(new_parameters);
    }
    
    return parameters;
}



#endif // GRADIENT_APPROXIMATION_HPP
