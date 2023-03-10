#include <vector>
#include <iostream>

#include "../inc/transfer_fcn.h"
#include "../inc/eigen_formatter.h"

namespace DT 
{
    
    namespace Examples
    {

        void example_transfer_function()
        {
            // input
            Eigen::VectorXd in = Eigen::VectorXd::Ones(500);

            // DC motor discrete fcn
            Eigen::VectorXd A {{ 1, -0.9802 }};
            Eigen::VectorXd B {{ 0.0594 }};
            DT::TransferFunction discrete_dc_model(B, A);
            discrete_dc_model.print();

            // simulation
            for (int i=0; i<499; i++)
            {
                double w = in[i];
                double y = discrete_dc_model.step(w);
                std::cout << y << std::endl;
            }
        } 
    }
}