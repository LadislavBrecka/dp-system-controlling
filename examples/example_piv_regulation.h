#include <vector>
#include <fstream>
#include <iostream>

#include "../inc/transfer_fcn.h"
#include "../inc/closed_loop.h"
#include "../inc/pole_placement.h"
#include "../inc/eigen_formatter.h"

namespace DT {

    namespace Examples
    {

        void example_piv_regulation() 
        {
            double T_step = 0.1;

            // input
            Eigen::VectorXd in = Eigen::VectorXd::Ones(500);

            // DC motor discrete fcn
            Eigen::VectorXd A {{ 1, -0.9802 }};
            Eigen::VectorXd B {{ 0.0594 }};
            DT::TransferFunction discrete_dc_model(B, A);
            discrete_dc_model.print();

            // convert discrete fcn to s-domain
            DT::TransferFunction continuous_dc_model;
            discrete_dc_model.d2c(T_step, continuous_dc_model);
            continuous_dc_model.print("s");

            // // make pole-placement with s-model
            auto PIV = DT::PolePlacement::PIV(continuous_dc_model, DT::TPZ, 2.0, 0.7, 1.0);  // omega = 2.0, b = 0.7, k = 1.0

            // // make PIV closed loop system with anti-windup algorithm
            DT::ClosedLoopSystem_PIV cls_piv(&discrete_dc_model, DT::TPZ, PIV.P, PIV.I, PIV.V, T_step, -2.0, 2.0, 1.0);

            std::cout << "Found PIV params: P: " << PIV.P << ", I: " << PIV.I << ", V: " << PIV.V << std::endl;

            // output to file for export to MATLAB
            std::ofstream y_log("logs/y_cpp.txt");
            std::ofstream u_log("logs/u_cpp.txt");
            std::ofstream e_log("logs/e_cpp.txt");
            y_log << "y_cpp = [\n"; u_log << "u_cpp = [\n"; e_log << "e_cpp = [\n";

            // P+IV regulator
            for (int i=0; i<499; i++)
            {
                double w = in[i];
                // double y = discrete_dc_model.step(w);
                DT::ClosedLoopStepResponse res = cls_piv.step(w);
                e_log << res.e << std::endl;
                u_log << res.u << std::endl;
                y_log << res.y << std::endl;
            }

            y_log << "]"; u_log << "]"; e_log << "]";
        }
    }
}