#pragma once

#include "./Eigen/Dense"
#include "./transfer_fcn.h"
#include "./common_types.h"

namespace DT
{
    namespace PolePlacement
    {
        // 0 zeros, 1 pole
        PIVRegCoefs PIV_0z_1p(const DT::TransferFunction& tf, DT::AproximationType aprox_type, double omega, double b, double k)
        {
            if (aprox_type == DT::PSD)
            {
#ifdef __x86_64__
                throw std::runtime_error("Cannot compute Pole Placement for PSD regulators!");
#else
                return {0,0,0};
#endif
            }

            const Eigen::VectorXd A = tf.get_denominator();
            const Eigen::VectorXd B = tf.get_numerator();

            // check conditions
            if (A.size() != 2 || B.size() != 1)
            {
#ifdef __x86_64__
                throw std::runtime_error("TF has wrong format!");
#else
                return {0,0,0};
#endif
            }

            double V = ((2*b*omega + k) * A[0] - 1) / B[0];
            double I = ((pow(omega, 2) + 2*b*omega*k) * A[0]) / B[0];
            double P = (pow(omega, 2)*k) / (pow(omega, 2) + 2*b*omega*k);

            return { P, I, V };
        }

        // 0 zeros, 2 poles
        PIDRegCoefs PID_0z_2p(const DT::TransferFunction& tf, DT::AproximationType aprox_type, double omega, double b, double k)
        {
            if (aprox_type == DT::PSD)
            {
#ifdef __x86_64__
                throw std::runtime_error("Cannot compute Pole Placement for PSD regulators!");
#else
                return {0,0,0};
#endif
            }

            const Eigen::VectorXd A = tf.get_denominator();
            const Eigen::VectorXd B = tf.get_numerator();

            // check conditions
            if (A.size() != 3 || B.size() != 1 || A[0] != 1.0)
            {
#ifdef __x86_64__
                throw std::runtime_error("TF has wrong format!");
#else
                return {0,0,0};
#endif
            }

            double P = (pow(omega, 2) + 2*b*omega*k - A[2]) / B[0];
            double I = (pow(omega, 2) * k) / B[0];
            double D = (2*b*omega + k - A[1]) / B[0];

            return { P, I, D };
        } 
    } 
}