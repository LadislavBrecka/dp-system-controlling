#pragma once

#include <cassert>
#include <iostream>

#include "./Eigen/Dense"
#include "./transfer_fcn.h"
#include "../inc/Exceptions/not_supported_exception.h"

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

    // TODO: need to implement POLE-PLACEMENT for PSD type regulators 
    namespace PolePlacement
    {

        PIVRegCoefs PIV(const DT::TransferFunction& tf, DT::AproximationType aprox_type, double omega, double b, double k)
        {
            if (aprox_type == DT::PSD)
                throw NotSupportedException( "PSD Aproximation type for pole-placement");

            const Eigen::VectorXd A = tf.get_denominator();
            const Eigen::VectorXd B = tf.get_numerator();

            // assert all conditions
            assert(A.size() == 2);

            double V = ((2*b*omega + k) * A[0] - 1) / B[0];
            double I = ((pow(omega, 2) + 2*b*omega*k) * A[0]) / B[0];
            double P = (pow(omega, 2)*k) / (pow(omega, 2) + 2*b*omega*k);

            return { P, I, V };
        }

        // TODO: need proper testing - not tested yet!!!
        PIDRegCoefs PID(const DT::TransferFunction& tf, DT::AproximationType aprox_type, double omega, double b, double k)
        {
            if (aprox_type == DT::PSD)
                throw NotSupportedException("PSD Aproximation type for pole-placement");

            const Eigen::VectorXd A = tf.get_denominator();
            const Eigen::VectorXd B = tf.get_numerator();

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