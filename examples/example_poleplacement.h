#include "../ControllingLib/inc/closed_loop.h"
#include "../ControllingLib/inc/transfer_fcn.h"
#include "../ControllingLib/inc/pole_placement.h"

namespace DT
{
    namespace Examples
    {
        void example_poleplacement()
        {
            double T_step = 0.1;

            // discrete fcn
            Eigen::VectorXd A{{1, -0.9802}};
            Eigen::VectorXd B{{0.0594}};
            DT::TransferFunction discrete_model(B, A);
            discrete_model.print();

            // convert discrete fcn to s-domain
            DT::TransferFunction continuous_model;
            discrete_model.d2c(T_step, continuous_model);
            continuous_model.print("s");

            // // make pole-placement with s-model
            auto PIV = DT::PolePlacement::PIV_0z_1p(continuous_model, DT::TPZ, 2.0, 0.7, 1.0); // omega = 2.0, b = 0.7, k = 1.0

            std::cout << "Found PIV params: P: " << PIV.P << ", I: " << PIV.I << ", V: " << PIV.V << std::endl;
        }
    }
}