#pragma once

#include "./Eigen/Dense"
#include <iostream>

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

    namespace PolePlacement
    {
        /*
        Works only for PID regulator with system of second order
        */
        RegCoefs evaluatePIDGains(const DT::TransferFunction& tf, double omega, double b, double k=0)
        {
            const Eigen::VectorXd A = tf.getDenominator();
            const Eigen::VectorXd B = tf.getNominator();
            const uint tf_order = A.size() - 1;

            // int P = 
            
            // if yes, make computations and return coeficients
            if (true)
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
}