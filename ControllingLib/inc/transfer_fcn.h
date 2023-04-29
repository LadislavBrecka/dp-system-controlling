#pragma once

#include <vector>
#include <memory>

#include "./Eigen/Dense"
#include "./Eigen/unsupported/Polynomials"
#include "./circle_buffer.h"

namespace DT 
{

    class TransferFunction;

    class TransferFunction
    {
    private:
        float ts;
        uint n_a, n_b;
        Eigen::VectorXd B, A;
        std::unique_ptr<DT::CircleBuffer> vU, vY;

    public:
        TransferFunction();
        TransferFunction(Eigen::VectorXd num, Eigen::VectorXd den);
        TransferFunction(std::vector<double> num, std::vector<double> den);
        ~TransferFunction();
        double step(double u);
        void d2c(double Ts,  DT::TransferFunction& c_tf);

#ifdef __x86_64__
        // NOT IMPLEMENTED IN STM32 PROJECT, IT'S JUST FOR DEBUGGING AND DEVELOPMENT
        void print(const std::string& var = "z");   
#endif

        inline Eigen::VectorXd get_numerator() const { return B; };
        inline Eigen::VectorXd get_denominator() const { return A; };

        void set_numerator(Eigen::VectorXd numerator);
        void set_numerator(std::vector<double> numerator);
        void set_denominator(Eigen::VectorXd denominator);
        void set_denominator(std::vector<double> denominator);
    };
}