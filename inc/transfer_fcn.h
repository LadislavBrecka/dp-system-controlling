#pragma once
#include "./Eigen/Dense"
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
        TransferFunction(Eigen::VectorXd nominator, Eigen::VectorXd denominator);
        TransferFunction(std::vector<double> nominator, std::vector<double> denominator);
        ~TransferFunction();
        double step(double u);
        void print();   
    };
    
}