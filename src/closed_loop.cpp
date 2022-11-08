#include "../inc/closed_loop.h"
#include "../inc/Eigen/Dense"
#include <iostream>

namespace DT 
{
    Regulator::Regulator()
    {
    }
    
    Regulator::~Regulator()
    {
    }

    void Regulator::initPsdRegulator(double p_gain, double i_gain, double d_gain)
    {
        Eigen::VectorXd i_nom {{ 1.0, 0.0 }};
        Eigen::VectorXd i_den {{ 1.0, -1.0 }};
        I_tf = std::make_unique<DT::TransferFunction>(i_nom, i_den); 

        Eigen::VectorXd d_den {{ 1.0, 0.0 }};
        Eigen::VectorXd d_nom {{ 1.0, -1.0 }};
        D_tf = std::make_unique<DT::TransferFunction>(d_nom, d_den); 

        P_gain = p_gain; I_gain = i_gain; D_gain = d_gain;
    }

    void Regulator::initDiscretePidRegulator(double N, double T, double p_gain, double i_gain, double d_gain)
    {
        Eigen::VectorXd i_nom {{ T/2.0, T/2.0 }};
        Eigen::VectorXd i_den {{ 1.0, -1.0 }};
        I_tf = std::make_unique<DT::TransferFunction>(i_nom, i_den); 

        Eigen::VectorXd d_nom {{ N, -N }};
        Eigen::VectorXd d_den {{ N*(T/2.0) + 1.0, N*(T/2.0) - 1.0 }};
        D_tf = std::make_unique<DT::TransferFunction>(d_nom, d_den); 

        P_gain = p_gain; I_gain = i_gain; D_gain = d_gain;
    }
    
    double Regulator::produceOutput(double e)
    {
        if (I_tf == nullptr || D_tf == nullptr) 
            throw std::domain_error("You must first call init() method for regulator to work properly!");
        
        double p_y = P_gain * e;
        double i_y = I_gain * I_tf->step(e);
        double d_y = D_gain * D_tf->step(e);
        return p_y + i_y + d_y;
    }
    
    ClosedLoopSystem::ClosedLoopSystem(DT::Regulator* reg, DT::TransferFunction* tf)
    {
        regulator = reg;
        system = tf;
    }
    
    ClosedLoopSystem::~ClosedLoopSystem()
    {
        regulator = nullptr;
        system = nullptr;
    }
    
    ClosedLoopStepResponse ClosedLoopSystem::step(double w)
    {
        double e = w - previous_y;
        double u = regulator->produceOutput(e);     

        double y = system->step(u);
        previous_y = y;
        return { e, u, y };
    }
}