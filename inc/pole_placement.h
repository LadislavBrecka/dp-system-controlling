#pragma once

namespace DT
{
    /*
    Need to implement pole placement method for PSD and discrete PID regulators (for now).

    */

    struct RegCoefs
    {
        double P;
        double I;
        double D;
    };
    
    RegCoefs evaluatePSD(double omega, double b, double k, const DT::TransferFunction& tf)
    {
        // check conditions, if system of linear equations has a solutions
        bool condition = true;

        // if yes, make computations and return coeficients
        if (condition)  
        {
            return { 1.0, 0.0, 0.0 };
        }
        // if no, throw exception
        else
        {
            throw std::domain_error("Regulator coeficients cannot be found by this method!");
        }
    }

    RegCoefs evaluateDiscretePID(double omega, double b, double k, const DT::TransferFunction& tf)
    {
        // check conditions, if system of linear equations has a solutions
        bool condition = true;
        
        // if yes, make computations and return coeficients
        if (condition)
        {
            return { 0.7, 0.8, 0.7 };
        }   
        // if no, throw exception
        else
        {
            throw std::domain_error("Regulator coeficients cannot be found by this method!");
        }
    }              
    
}