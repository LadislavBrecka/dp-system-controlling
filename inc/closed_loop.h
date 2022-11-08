#pragma once

#include "./transfer_fcn.h"
#include <memory>

/*---------- Example of usage in main.cpp: ------------

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

    DT::Regulator reg;
    reg.initDiscretePidRegulator(10, 0.1, 0.7, 0.8, 0.7);

    DT::ClosedLoopSystem cls(&reg, &tf);

    for (int i=0; i<500; i++)
    {
        double w = in[i];
        DT::ClosedLoopStepResponse out = cls.step(w);
        std::cout << out.y << std::endl;
    }

    return 0;
}
 
 
*/

namespace DT
{
    struct ClosedLoopStepResponse
    {
        double e;
        double u;
        double y;
    };

    class Regulator 
    {
    private:
        double P_gain;
        double I_gain;
        double D_gain;
        std::unique_ptr<DT::TransferFunction> I_tf;
        std::unique_ptr<DT::TransferFunction> D_tf;

    public:
        Regulator();
        ~Regulator();
    
        void initPsdRegulator(double p_gain, double i_gain, double d_gain);
        void initDiscretePidRegulator(double N, double T, double p_gain, double i_gain, double d_gain);
        double produceOutput(double e);       
    };

    class ClosedLoopSystem
    {
    private:
        DT::Regulator* regulator;
        DT::TransferFunction* system; 
        double previous_y = 0;

    public:
        ClosedLoopSystem(DT::Regulator* reg, DT::TransferFunction* tf);
        ~ClosedLoopSystem();
        ClosedLoopStepResponse step(double w);
    };
}
