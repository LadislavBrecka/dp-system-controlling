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
    ClosedLoopSystem_PIV::ClosedLoopSystem_PIV(DT::TransferFunction* tf, AproximationType aproxType, double P, double I, double V, double T)
    : P_gain(P), I_gain(I), V_gain(V)
    {
        position_integrator = std::make_unique<Integrator>(aproxType, T);
        i_reg_integrator = std::make_unique<Integrator>(aproxType, T);
        system = tf;
    }
    
    ClosedLoopSystem_PIV::~ClosedLoopSystem_PIV()
    {
        system = nullptr;
    }
    
    ClosedLoopStepResponse ClosedLoopSystem_PIV::step(double w)
    {
        double e_pos = w - previous_iy;
        double u_pos = P_gain * e_pos;
    
        double e_speed = u_pos - previous_y;
        double u_speed = i_reg_integrator->produceOutput(e_speed) * I_gain;

        double u = u_speed - previous_y * V_gain;

        double speed = system->step(u);
        double position = position_integrator->produceOutput(speed);
        previous_y = speed;
        previous_iy = position;

        return { e_pos, u, position };
    }
}