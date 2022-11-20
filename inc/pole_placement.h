#pragma once

#include "./Eigen/Dense"
#include "./transfer_fcn.h"
#include <iostream>
#include <cassert>

namespace DT
{
    struct PIDRegCoefs
    {
        double P;
        double I;
        double D;
    };

    struct PIVRegCoefs
    {
        double P;
        double I;
        double V;
    };

    namespace PolePlacement
    {
        PIVRegCoefs PIV(const DT::TransferFunction& tf, double omega, double b, double k)
        {
            const Eigen::VectorXd A = tf.getDenominator();
            const Eigen::VectorXd B = tf.getNominator();

            // assert all conditions
            assert(A.size() == 2);

            double V = ((2*b*omega + k) * A[0] - 1) / B[0];
            double I = ((pow(omega, 2) + 2*b*omega*k) * A[0]) / B[0];
            double P = (pow(omega, 2)*k) / (pow(omega, 2) + 2*b*omega*k);

            return { P, I, V };
        }

        // TODO: make pole placement algorithm for any PID regulators
        // PIDRegCoefs PID(const DT::TransferFunction& tf, double beta, double gama, double delta)
        // {
        //     const Eigen::VectorXd A = tf.getDenominator();
        //     const Eigen::VectorXd B = tf.getNominator();
        //     const uint n_a = A.size() - 1;

        //     assert(B.size() == 1);
        //     assert(A.size() == 5);

        //     int P = gama / B[0];
        //     int I = delta / B[0];
        //     int D = (beta - A[1]) / B[0];
            
        //     // if yes, make computations and return coeficients
        //     if (true)
        //     {
        //         return { P, I, D };
        //     }   
        //     // if no, throw exception
        //     else
        //     {
        //         throw std::domain_error("PIDRegulator coeficients cannot be found by this method!");
        //     }
        // } 
    } 
}