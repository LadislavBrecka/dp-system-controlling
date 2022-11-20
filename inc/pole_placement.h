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

        // TODO: need proper testing - not tested yet!!!
        PIDRegCoefs PID(const DT::TransferFunction& tf, double omega, double b, double k)
        {
            const Eigen::VectorXd A = tf.getDenominator();
            const Eigen::VectorXd B = tf.getNominator();

            assert(B.size() == 1);
            assert(A.size() == 3);
            assert(A[0] == 1.0);

            double P = (pow(omega, 2) + 2*b*omega*k - A[2]) / B[0];
            double I = (pow(omega, 2) * k) / B[0];
            double D = (2*b*omega + k - A[1]) / B[0];

            return { P, I, D };
        } 
    } 
}