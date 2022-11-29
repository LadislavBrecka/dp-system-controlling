#include "inc/identification.h"
#include "inc/transfer_fcn.h"
#include "inc/closed_loop.h"
#include "inc/pole_placement.h"
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm> 

const Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");

int main(int argc, char const *argv[])
{
    double T_step = 0.1;

    // input
    Eigen::VectorXd in = Eigen::VectorXd::Ones(500);

    // DC motor discrete fcn
    Eigen::VectorXd A {{ 1.0, -0.9802 }};
    Eigen::VectorXd B {{ 0.0594 }};
    DT::TransferFunction discrete_dc_model(B, A);

    Eigen::VectorXd cA {{ 5, 1 }};
    Eigen::VectorXd cB {{ 3 }};
    DT::TransferFunction continuous_dc_model(cB, cA);

    // make pole-placement 
    auto PIV = DT::PolePlacement::PIV(continuous_dc_model, DT::PSD, 2.0, 0.7, 1.0);  // omega = 2.0, b = 0.7, k = 1.0

    // make PIV closed loop system with anti-windup algorithm
    DT::ClosedLoopSystem_PIV cls_piv(&discrete_dc_model, DT::TPZ, PIV.P, PIV.I, PIV.V, T_step, -0.1, 0.1, 1.0);

    std::cout << "Found PIV params: P: " << PIV.P << ", I: " << PIV.I << ", V: " << PIV.V << std::endl;

    // P+IV regulator
    for (int i=0; i<499; i++)
    {
        double w = in[i];
        // double y = discrete_dc_model.step(w);
        DT::ClosedLoopStepResponse res = cls_piv.step(w);
        std::cout << res.y << std::endl; 
    }

    return 0;
}