#pragma once

#include "./transfer_fcn.h"
#include <memory>

namespace DT
{
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
    
        void init_psd(double p_gain, double i_gain, double d_gain);
        void init_discrete_pid(double N, double T, double p_gain, double i_gain, double d_gain);
        double produce_output(double e);       
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
        double step(double w);
    };
}
