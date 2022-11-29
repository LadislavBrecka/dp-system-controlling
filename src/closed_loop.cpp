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
    PIDRegulator::PIDRegulator(AproximationType aproxType, double P, double I, double D, double T, double N,
                               double uMin, double uMax, double Kaw)
    : P_gain(P), I_gain(I), D_gain(D), u_min(uMin), u_max(uMax), k_aw(Kaw)
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
        double i_y = integrator->produceOutput(e * I_gain + prev_aw_gain);
        double d_y = derivator->produceOutput(e * D_gain);

        double u = p_y + i_y + d_y;

        // anti-wind up algorithm
        double u_before_saturation = u;
        if (u > u_max)          u = u_max;
        else if (u < u_min)     u = u_min;
        double diff = u - u_before_saturation;
        double aw_gain = diff * k_aw;
        prev_aw_gain = aw_gain;  

        return u;
    }

    /*
    PID SYSTEM IMPLEMENTATIONS
    */
    ClosedLoopSystem_PID::ClosedLoopSystem_PID(DT::TransferFunction* tf, AproximationType aproxType, 
                                                double P, double I, double D, double T, double N,
                                                double uMin, double uMax, double Kaw)
    {
        pid_regulator = std::make_unique<DT::PIDRegulator>(aproxType, P, I, D, T, N, uMin, uMax, Kaw);
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
    ClosedLoopSystem_PIV::ClosedLoopSystem_PIV(DT::TransferFunction* tf, AproximationType aproxType, 
                                                double P, double I, double V, double T,
                                                double uMin, double uMax, double Kaw)
    : P_gain(P), I_gain(I), V_gain(V), u_min(uMin), u_max(uMax), k_aw(Kaw)
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
        double e_1 = w - previous_iy;                                                   // position error
        double u_1 = P_gain * e_1;                                                      // position correction signal
    
        double e_2 = u_1 - previous_y;                                                  // speed error
        double u_2 = i_reg_integrator->produceOutput(e_2 * I_gain + prev_aw_gain);      // speed correction signal

        double u = u_2 - previous_y * V_gain;                                           // final correction signal

        // anti-wind up algorithm
        double u_before_saturation = u;
        if (u > u_max)          u = u_max;
        else if (u < u_min)     u = u_min;
        double diff = u - u_before_saturation;
        double aw_gain = diff * k_aw;
        prev_aw_gain = aw_gain;  

        double y = system->step(u);                                                     // y  - speed
        double iy = output_integrator->produceOutput(y);                                // iy - position
        previous_y = y;
        previous_iy = iy;

        return { e_1, u, iy };
    }
}