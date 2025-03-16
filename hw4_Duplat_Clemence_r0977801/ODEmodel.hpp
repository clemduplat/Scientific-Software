#ifndef ODEMODEL_HPP
#define ODEMODEL_HPP

#include <vector>
#include <functional>
#include <algorithm>  
#include <iterator>
#include <numeric>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

class ODEModel {
public:
    ODEModel(size_t systemSize)
        : systemSize(systemSize) {}

    vector<double> calculateSystemDerivatives(const vector<double>& x) const {
        vector<double> derivatives(systemSize);
        size_t idx = 0; // Initialize index
        std::transform(x.begin(), x.end(), derivatives.begin(), 
            [this, &idx](double xi) mutable {
                double result = -10 * std::pow(xi - 0.1 * idx, 3);
                ++idx; // Increment index
                return result;
            });
        return derivatives;
    }

    matrix<double> calculateJacobian(const vector<double>& x) const {
        matrix<double> jacobian(systemSize, systemSize, 0.0); 

        std::vector<size_t> indices(systemSize);
        std::iota(indices.begin(), indices.end(), 0);

        std::transform(indices.begin(), indices.end(), indices.begin(),
            [&jacobian, &x](size_t i) {
                jacobian(i, i) = -30 * std::pow(x(i) - 0.1 * i, 2);
                return i; 
            });

        return jacobian;
    }

private:
    size_t systemSize;
};

#endif // ODEMODEL_HPP