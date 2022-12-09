#pragma once
#include "./Eigen/Dense"
// #include "./unsupported/Eigen/Polynomials"
#include <vector>
#include <memory>

#include "circle_buffer.h"

namespace DT {

    class TransferFunction;

    class TransferFunction
    {
    private:
        float ts;
        Eigen::VectorXd B;
        Eigen::VectorXd A;
        uint n_a; uint n_b;

        std::unique_ptr<DT::CircleBuffer> vU;
        std::unique_ptr<DT::CircleBuffer> vY;

    public:
        TransferFunction();
        TransferFunction(Eigen::VectorXd nominator, Eigen::VectorXd denominator);
        TransferFunction(std::vector<double> nominator, std::vector<double> denominator);
        ~TransferFunction();
        double step(double u);
        void d2c(double Ts,  DT::TransferFunction& c_tf);
        void print(const std::string& var = "z");   

        inline uint getPolynomialOrder() const { return n_a; };
        inline Eigen::VectorXd getNominator() const { return B; };
        inline Eigen::VectorXd getDenominator() const { return A; };

        void setNominator(Eigen::VectorXd nominator);
        void setNominator(std::vector<double> nominator);
        void setDenominator(Eigen::VectorXd denominator);
        void setDenominator(std::vector<double> denominator);
    };
}