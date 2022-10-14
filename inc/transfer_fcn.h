#pragma once
#include "./Eigen/Dense"
#include <vector>

namespace DT {

    class TransferFunction;

    class TransferFunction
    {
    private:
        float ts;
        Eigen::VectorXd B;
        Eigen::VectorXd A;
        int n_a; int n_b;  // need to be int

        Eigen::VectorXd vU;
        Eigen::VectorXd vY;

    public:
        TransferFunction(Eigen::VectorXd nominator, Eigen::VectorXd denominator);
        TransferFunction(std::vector<double> nominator, std::vector<double> denominator);
        ~TransferFunction();
        double step(double u);
        void print();

    private:
        void shiftU(double u);
        void shiftY(double y);
    };
   
}