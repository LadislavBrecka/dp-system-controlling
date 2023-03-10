#pragma once

#include <memory>

#include "./transfer_fcn.h"

// structures and custom data types
namespace DT
{
    struct ClosedLoopStepResponse
    {
        double e;
        double u;
        double y;
    };

    enum AproximationType
    {
        FWD = 0,
        BWD = 1,
        TPZ = 2,
        PSD = 3
    };
}

// blocks for regulators
namespace DT
{
    class Integrator : public TransferFunction
    {
    public:
        Integrator(AproximationType aprox_type, double T=0.0);
        ~Integrator();
    };

    class Derivator : public TransferFunction
    {
    public:
        Derivator(AproximationType aprox_type, double T=0.0, double N=0.0);
        ~Derivator();
    };

    class PIDRegulator 
    {
    private:
        double P_gain, I_gain, D_gain;
        double u_min, u_max, k_aw, prev_aw_gain = 0.0;
        std::unique_ptr<Integrator> integrator;
        std::unique_ptr<Derivator> derivator;

    public:
        PIDRegulator(AproximationType aprox_type, double P, double I, double D, double T=0.0, double N=0.0,
                     double uMin=-10.0, double uMax=10.0, double Kaw=0.0);
        ~PIDRegulator();

        double step(double e);       
    };
}

// specifying closed loop systems
namespace DT 
{

    class ClosedLoopSystem_PID
    {
    private:
        std::unique_ptr<DT::PIDRegulator> pid_regulator;
        DT::TransferFunction* system; 
        double previous_y = 0.0;

    public:
        ClosedLoopSystem_PID(DT::TransferFunction* tf, AproximationType aprox_type, 
                             double P, double I, double D, double T=0.0, double N=0.0,
                             double uMin=-10.0, double uMax=10.0, double Kaw=0.0);
        ~ClosedLoopSystem_PID();
        ClosedLoopStepResponse step(double w);
    };

    class ClosedLoopSystem_PIV
    {
    private:
        double P_gain, I_gain, V_gain;
        double u_min, u_max, k_aw, prev_aw_gain = 0.0;
        std::unique_ptr<DT::Integrator> output_integrator;
        std::unique_ptr<DT::Integrator> i_reg_integrator;
        DT::TransferFunction* system; 
        double previous_y = 0.0;
        double previous_iy = 0.0;

    public:
        ClosedLoopSystem_PIV(DT::TransferFunction* tf, AproximationType aprox_type, 
                             double P, double I, double V, double T=0.0,
                             double uMin=-10.0, double uMax=10.0, double Kaw=0.0);
        ~ClosedLoopSystem_PIV();
        ClosedLoopStepResponse step(double w);
    };
}