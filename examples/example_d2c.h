#include "../inc/transfer_fcn.h"

namespace DT
{
    namespace Examples
    {
        void example_d2c()
        {
            double T_step = 0.1;

            // discrete fcn
            Eigen::VectorXd A{{1, 0.6, 1.2}};
            Eigen::VectorXd B{{2}};
            DT::TransferFunction discrete_model(B, A);
            discrete_model.print();

            // convert discrete fcn to s-domain
            DT::TransferFunction continuous_model;
            discrete_model.d2c(T_step, continuous_model);
            continuous_model.print("s");
        }
    }
}