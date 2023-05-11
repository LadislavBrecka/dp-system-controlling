#pragma once

#include <memory>

#include "./Eigen/Dense"
#include "./common_types.h"
#include "./transfer_fcn.h"

namespace DT
{
    class Integrator : public DT::TransferFunction
    {
    public:
        Integrator(DT::AproximationType aprox_type, double T);
        ~Integrator();
    };

    class Derivator : public DT::TransferFunction
    {
    public:
        Derivator(DT::AproximationType aprox_type, double T, double N);
        ~Derivator();
    };
}

// declarations of different regulators
namespace DT
{
    class PIDRegulator 
    {
    private:
        double P_gain, I_gain, D_gain;
        double u_min, u_max, k_aw, prev_aw_gain = 0.0;
        std::unique_ptr<DT::Integrator> integrator;
        std::unique_ptr<DT::Derivator> derivator;

    public:
        PIDRegulator(DT::AproximationType aprox_type, double P, double I, double D, double T, double N,
                     double uMin=-10.0, double uMax=10.0, double Kaw=0.0);
        ~PIDRegulator();

        DT::RegulatorResponse step(double w, double prev_y);      
    };

    class PIVRegulator
    {
    private:
        double P_gain, I_gain, V_gain;
        double u_min, u_max, k_aw, prev_aw_gain = 0.0;
        std::unique_ptr<DT::Integrator> i_reg_integrator;

    public:
        PIVRegulator(DT::AproximationType aprox_type, double P, double I, double V, double T,
                     double uMin=-10.0, double uMax=10.0, double Kaw=0.0);
        ~PIVRegulator();

        DT::RegulatorResponse step(double w, double prev_y, double prev_iy);
    };
}

// declarations of different closed loop systems
namespace DT 
{
    class ClosedLoopSystem_PID
    {
    private:
        std::unique_ptr<DT::PIDRegulator> pid_regulator;
        DT::TransferFunction* system; 
        double prev_y = 0.0;

    public:
        ClosedLoopSystem_PID(DT::TransferFunction* tf, AproximationType aprox_type, 
                             double P, double I, double D, double T, double N,
                             double uMin=-10.0, double uMax=10.0, double Kaw=0.0);
        ~ClosedLoopSystem_PID();
        DT::ClosedLoopStepResponse step(double w);
    };

    class ClosedLoopSystem_PIV
    {
    private:
        std::unique_ptr<DT::Integrator> output_integrator;
        std::unique_ptr<DT::PIVRegulator> piv_regulator;
        DT::TransferFunction* system; 
        double prev_y = 0.0;
        double prev_iy = 0.0;

    public:
        ClosedLoopSystem_PIV(DT::TransferFunction* tf, DT::AproximationType aprox_type, 
                             double P, double I, double V, double T,
                             double uMin=-10.0, double uMax=10.0, double Kaw=0.0);
        ~ClosedLoopSystem_PIV();
        DT::ClosedLoopStepResponse step(double w);
    };
}