#include "../inc/transfer_fcn.h"
#include <iostream>
#include <cmath>

namespace DT {

    TransferFunction::TransferFunction(Eigen::VectorXd nominator, Eigen::VectorXd denominator)
    {
        n_a = denominator.size();
        n_b = nominator.size();

        A = denominator;
        B = nominator;

        vU = std::make_unique<DT::CircleBuffer>(Eigen::VectorXd::Zero(n_b));
        vY = std::make_unique<DT::CircleBuffer>(Eigen::VectorXd::Zero(n_a));
    }
    
    TransferFunction::TransferFunction(std::vector<double> nominator, std::vector<double> denominator)
    {
        n_a = denominator.size();
        n_b = nominator.size();

        A = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(denominator.data(), n_a);
        B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(nominator.data(), n_b);

        vU = std::make_unique<DT::CircleBuffer>(Eigen::VectorXd::Zero(n_b));
        vY = std::make_unique<DT::CircleBuffer>(Eigen::VectorXd::Zero(n_a));
    }
    
    TransferFunction::~TransferFunction()
    {
    }
    
    double TransferFunction::step(double u)
    {
        // shift vector u so it contains the newest input sample
        vU->add(u);

        // input part
        double input_part = 0;
        for (uint i = 0; i < n_b; i++)
        {
            // input_part += B[i] * vU->operator[](i);
            input_part += B[i] * vU->at(i);
        }

        double output_part = 0;
        for (uint i = 1; i < n_a; i++)
        {
            output_part += A[i] * vY->at(i-1);
        }

        double sum = input_part + output_part;
        double y = 0.0;
        if (sum != 0.0) {
            y = (1 / A[0]) * (input_part - output_part);
        }
        
        // shift vector y so it containts the newest output sample
        vY->add(y);
        return y;
    }
    
    void TransferFunction::print()
    {
        // nominator printing
        for (uint i = 0; i < n_b; i++)
        {
            if (B[i] != 0.0)
            {
                std::string z_index = (i == 0 ? "" : "z-" + std::to_string(i));
                std::cout << fabs(B[i]) << z_index;
                if (i != n_b - 1)
                {
                    if (B[i] > 0.0) std::cout << " + "; else std::cout << " - ";
                }    
            }
        }

        // dividing line printing
        std::cout << std::endl;
        for (uint i = 0; i < std::max(n_a, n_b) * 8; i++)
        {
            std::cout << "-";
        }
        std::cout << std::endl;

        // denominator printing
        for (uint i = 0; i < n_a; i++)
        {
            if (A[i] != 0.0)
            {
                std::string z_index = (i == 0 ? "" : "z-" + std::to_string(i));
                std::cout << fabs(A[i]) << z_index;
                if (i != n_a - 1)
                {
                    if (A[i] > 0.0) std::cout << " + "; else std::cout << " - ";
                }
            }
        }

        // new line at the end
        std::cout << std::endl;
    }

}