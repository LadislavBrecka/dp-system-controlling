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
    std::ifstream is("input_signal.txt");
    std::istream_iterator<double> start(is), end;
    std::vector<double> in(start, end);
    std::cout << "Read " << in.size() << " numbers" << std::endl;

    Eigen::Vector3d A {{ 1.0, -1.8955, 0.9050 }};
    Eigen::Vector3d B {{ 0.0024, 0.0048, 0.0024 }};
    DT::TransferFunction tf(B, A);

    // there will be poleplacement method
    DT::RegCoefs coefs = DT::PolePlacement::evaluatePIDGains(tf, 2, 0.7, 2.0);

    DT::Regulator reg;
    reg.initDiscretePidRegulator(10, 0.1, coefs.P, coefs.I, coefs.D);

    DT::ClosedLoopSystem cls(&reg, &tf);

    // for (int i=0; i<500; i++)
    // {
    //     double w = in[i];
    //     DT::ClosedLoopStepResponse out = cls.step(w);
    //     std::cout << out.y << std::endl;
    // }

    return 0;
}
