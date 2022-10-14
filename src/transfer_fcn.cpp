#include "../inc/transfer_fcn.h"
#include <iostream>

namespace DT {

    TransferFunction::TransferFunction(Eigen::VectorXd nominator, Eigen::VectorXd denominator)
    {
        n_a = denominator.size();
        n_b = nominator.size();

        A = denominator;
        B = nominator;

        vU = Eigen::VectorXd::Zero(n_b);
        vY = Eigen::VectorXd::Zero(n_a);
    }
    
    TransferFunction::TransferFunction(std::vector<double> nominator, std::vector<double> denominator)
    {
        n_a = denominator.size();
        n_b = nominator.size();

        A = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(denominator.data(), n_a);
        B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(nominator.data(), n_b);

        vU = Eigen::VectorXd::Zero(n_b);
        vY = Eigen::VectorXd::Zero(n_a);
    }
    
    TransferFunction::~TransferFunction()
    {
    }
    
    double TransferFunction::step(double u)
    {
        // shift vector u so it contains the newest input sample
        shiftU(u);

        // input part
        double input_part = 0;
        for (size_t i = 0; i < n_b; i++)
        {
            input_part += B[i] * vU[i];
        }

        double output_part = 0;
        for (size_t i = 1; i < n_a; i++)
        {
            output_part += A[i] * vY[i-1];
        }

        double sum = input_part + output_part;
        double y = 0.0;
        if (sum != 0.0) {
            y = (1 / A[0]) * (input_part - output_part);
        }
        
        // shift vector y so it containts the newest output sample
        shiftY(y);
        return y;
    }
    
    void TransferFunction::print()
    {
        std::cout << "Not implemented yet!" << std::endl;
    }
    
    void TransferFunction::shiftU(double u)
    {
        // update vector u with inputs of system
        for (size_t i = n_b-1; i >= 1; i--) vU[i] = vU[i-1]; 
        vU[0] = u;
    }
    
    void TransferFunction::shiftY(double y)
    {
         // update vector u with inputs of system
        for (size_t i = n_a-1; i >= 1; i--) vY[i] = vY[i-1]; 
        vY[0] = y;
    }


}