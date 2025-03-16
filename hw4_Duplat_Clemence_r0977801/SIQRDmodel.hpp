#ifndef SIQRDMODEL_HPP
#define SIQRDMODEL_HPP
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <functional>
using namespace boost::numeric::ublas;
struct SIQRDParams {
    double beta;
    double mu;
    double gamma;
    double alpha;
    double delta;
};

class SIQRDModel{
    public: 

        SIQRDModel(SIQRDParams params, vector<double> initial_state)
        : params(params), state(initial_state) {}

        void set_params(const SIQRDParams& new_params) {
            params = new_params;
        }

        SIQRDParams get_params() const {
            return params;
        }

        void set_initial_state(const vector<double>& new_state) {
            state = new_state;
        }

        vector<double> get_initial_state() const {
            return state;
        }

        vector<double> calculate_derivatives(const vector<double>& x) const{
            double S=x(0);
            double I=x(1);
            double Q=x(2);
            double R=x(3);
            //double D=x(4);
            vector<double> derivatives(5);
            derivatives(0) = S_der(S, I, R);
            derivatives(1) = I_der(S, I, R);
            derivatives(2) = Q_der(I, Q);
            derivatives(3) = R_der(I, Q, R);
            derivatives(4) = D_der(I, Q);
            return derivatives;
        }
        matrix<double> calculate_jacobian(const vector<double>& x) const {
            double S=x(0); 
            double I=x(1); 
            double R=x(3);
            matrix<double> Jacobian(5, 5, 0.0);  

            Jacobian(0,0) = (params.beta * I /(S+I+R))*((S/(S+I+R)) -1 );
            Jacobian(0,1) = (params.beta * S /(S+I+R) )*((I/(S+I+R)) -1);
            Jacobian(0,3) = params.beta * S *I/(std::pow((S + I + R),2)) + params.mu;

            Jacobian(1,0) = (params.beta * I/(S+I+R)) * (1-S/(S+I+R));
            Jacobian(1,1) = (params.beta*S /(S+I+R))*(1-I/(S+I+R)) - params.gamma - params.delta - params.alpha;
            Jacobian(1,3)= - (params.beta * I* S) / (std::pow((S + I + R),2));

            Jacobian(2,1) = params.delta;
            Jacobian(2,2) = -(params.gamma + params.alpha);

            Jacobian(3,1) = params.gamma;
            Jacobian(3,2) = params.gamma;
            Jacobian(3,3) = -params.mu;

            Jacobian(4,1) = params.alpha;
            Jacobian(4,2) = params.alpha;

            return Jacobian;
        }
    private:
        SIQRDParams params;
        vector<double> state;
        double S_der(double S, double I, double R) const {
            return -params.beta * (I / (S + I + R)) * S + params.mu * R;
        }

        double I_der(double S, double I, double R) const{
            return (params.beta * (S / (S + I + R)) - params.gamma - params.delta - params.alpha) * I;
        }

        double Q_der(double I, double Q) const{
            return params.delta * I - (params.gamma + params.alpha) * Q;
        }

        double R_der(double I, double Q, double R) const{
            return params.gamma * (I + Q) - params.mu * R;
        }

        double D_der(double I, double Q) const{
            return params.alpha * (I + Q);
        }

};
#endif