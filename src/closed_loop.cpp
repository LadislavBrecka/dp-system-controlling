#include "../inc/closed_loop.h"
#include "../inc/Eigen/Dense"
#include <iostream>

namespace DT 
{

    /*
    DERIVATOR AND INTEGRATOR IMPLEMENTATIONS
    */
   Integrator::Integrator(AproximationType type, double T)
   {
        switch(type)
        {
            case TPZ:
            {
                Eigen::VectorXd i_nom {{ T, T }};
                Eigen::VectorXd i_den {{ 2.0, -2.0 }};
                tf = std::make_unique<DT::TransferFunction>(i_nom, i_den); 
                break;
            }

            case PSD:
            {
                Eigen::VectorXd i_nom {{ 1.0, 0.0 }};
                Eigen::VectorXd i_den {{ 1.0, -1.0 }};
                tf = std::make_unique<DT::TransferFunction>(i_nom, i_den);
                break;
            }
        }
   }

   Integrator::~Integrator()
   {
   }

   Derivator::Derivator(AproximationType type, double T, double N)
   {
        switch(type)
        {
            case TPZ:
            {
                Eigen::VectorXd d_nom {{ N, -N }};
                Eigen::VectorXd d_den {{ N*(T/2.0) + 1.0, N*(T/2.0) - 1.0 }};
                tf = std::make_unique<DT::TransferFunction>(d_nom, d_den); 
                break;
            }

            case PSD:
            {
                Eigen::VectorXd d_den {{ 1.0, 0.0 }};
                Eigen::VectorXd d_nom {{ 1.0, -1.0 }};
                tf = std::make_unique<DT::TransferFunction>(d_nom, d_den);
                break;
            }
        }
   }

   Derivator::~Derivator()
   {
   }

    /*
    REGULATOR IMPLEMENTATIONS
    */
    PIDRegulator::PIDRegulator(AproximationType aproxType, double P, double I, double D, double T, double N)
    : P_gain(P), I_gain(I), D_gain(D)
    {
        integrator = std::make_unique<Integrator>(aproxType, T);
        derivator = std::make_unique<Derivator>(aproxType, N, T);
    }
    
    PIDRegulator::~PIDRegulator()
    {
    }
    
    double PIDRegulator::produceOutput(double e)
    {
        if (integrator == nullptr || derivator == nullptr) 
            throw std::domain_error("You must first call init() method for regulator to work properly!");
        
        double p_y = P_gain * e;
        double i_y = I_gain * integrator->produceOutput(e);
        double d_y = D_gain * derivator->produceOutput(e);
        return p_y + i_y + d_y;
    }

    /*
    PID SYSTEM IMPLEMENTATIONS
    */
    ClosedLoopSystem_PID::ClosedLoopSystem_PID(DT::TransferFunction* tf, AproximationType aproxType, double P, double I, double D, double T, double N)
    {
        pid_regulator = std::make_unique<DT::PIDRegulator>(aproxType, P, I, D, T, N);
        system = tf;
    }
    
    ClosedLoopSystem_PID::~ClosedLoopSystem_PID()
    {
        system = nullptr;
    }
    
    ClosedLoopStepResponse ClosedLoopSystem_PID::step(double w)
    {
        double e = w - previous_y;
        double u = pid_regulator->produceOutput(e);     

        double y = system->step(u);
        previous_y = y;
        return { e, u, y };
    }

    /*
    PIV SYSTEM IMPLEMENTATIONS
    */
    ClosedLoopSystem_PIV::ClosedLoopSystem_PIV(DT::TransferFunction* tf, AproximationType aproxType, double P, double I, double V, double T=0.0)
    : P_gain(P), I_gain(I), V_gain(V)
    {
        output_integrator = std::make_unique<Integrator>(aproxType, T);
        i_reg_integrator = std::make_unique<Integrator>(aproxType, T);
        system = tf;
    }
    
    ClosedLoopSystem_PIV::~ClosedLoopSystem_PIV()
    {
        system = nullptr;
    }
    
    ClosedLoopStepResponse ClosedLoopSystem_PIV::step(double w)
    {
        double e_1 = w - previous_iy;                       // position error
        double u_1 = P_gain * e_1;                          // position correction signal
    
        double e_2 = u_1 - previous_y;                                  // speed error
        double u_2 = i_reg_integrator->produceOutput(e_2) * I_gain;     // speed correction signal

        double u = u_2 - previous_y * V_gain;               // final correction signal

        double y = system->step(u);                         // y  - speed
        double iy = output_integrator->produceOutput(y);    // iy - position
        previous_y = y;
        previous_iy = iy;

        return { e_1, u, iy };
    }
}