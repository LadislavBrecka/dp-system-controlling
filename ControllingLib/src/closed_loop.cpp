#include "../inc/closed_loop.h"

namespace DT
{
    /*
    INTEGRATOR IMPLEMENTATIONS
    */
    Integrator::Integrator(DT::AproximationType aprox_type, double T)
    {
        switch (aprox_type)
        {
            default:
            case DT::TPZ:
                set_numerator(Eigen::VectorXd{{T, T}});
                set_denominator(Eigen::VectorXd{{2.0, -2.0}});
                break;

            case DT::PSD:
                set_numerator(Eigen::VectorXd{{1.0, 0.0}});
                set_denominator(Eigen::VectorXd{{1.0, -1.0}});
                break;
        }
    }

    Integrator::~Integrator()
    {
    }

    /*
    DERIVATOR IMPLEMENTATIONS
    */
    Derivator::Derivator(DT::AproximationType aprox_type, double T, double N)
    {
        switch (aprox_type)
        {
            default:
            case DT::TPZ:
                set_numerator(Eigen::VectorXd{{N, -N}});
                set_denominator(Eigen::VectorXd{{N * (T / 2.0) + 1.0, N * (T / 2.0) - 1.0}});
                break;

            case DT::PSD:
                set_numerator(Eigen::VectorXd{{1.0, -1.0}});
                set_denominator(Eigen::VectorXd{{1.0, 0.0}});
                break;
        }
    }

    Derivator::~Derivator()
    {
    }
}


namespace DT
{
    /*
    PID REGULATOR IMPLEMENTATIONS
    */
    PIDRegulator::PIDRegulator(DT::AproximationType aproxType, double P, double I, double D, double T, double N,
                               double uMin, double uMax, double Kaw)
        : P_gain(P), I_gain(I), D_gain(D), u_min(uMin), u_max(uMax), k_aw(Kaw)
    {
        integrator = std::make_unique<DT::Integrator>(aproxType, T);
        derivator = std::make_unique<DT::Derivator>(aproxType, N, T);
    }

    PIDRegulator::~PIDRegulator()
    {
    }

    DT::RegulatorResponse PIDRegulator::step(double w, double prev_y)
    {
        double e = w - prev_y;
    
        double p_y = P_gain * e;
        double i_y = integrator->step(e * I_gain + prev_aw_gain);
        double d_y = derivator->step(e * D_gain);

        double u = p_y + i_y + d_y;

        // anti-wind up algorithm
        double u_before_saturation = u;
        if (u > u_max)
            u = u_max;
        else if (u < u_min)
            u = u_min;
        double diff = u - u_before_saturation;
        double aw_gain = diff * k_aw;
        prev_aw_gain = aw_gain;

        return { u, e };
    }

     /*
    PIV REGULATOR IMPLEMENTATIONS
    */
    PIVRegulator::PIVRegulator(DT::AproximationType aproxType, double P, double I, double V, double T,
                               double uMin, double uMax, double Kaw)
        : P_gain(P), I_gain(I), V_gain(V), u_min(uMin), u_max(uMax), k_aw(Kaw)
    {
        i_reg_integrator = std::make_unique<DT::Integrator>(aproxType, T);
    }

    PIVRegulator::~PIVRegulator()
    {
    }

    DT::RegulatorResponse PIVRegulator::step(double w, double prev_y, double prev_iy)
    {
        double e_1 = w - prev_iy; // position error
        double u_1 = P_gain * e_1;    // position correction signal

        double e_2 = u_1 - prev_y;                                    // speed error
        double u_2 = i_reg_integrator->step(e_2 * I_gain + prev_aw_gain); // speed correction signal

        double u = u_2 - prev_y * V_gain; // final correction signal

        // anti-wind up algorithm
        double u_before_saturation = u;
        if (u > u_max)
            u = u_max;
        else if (u < u_min)
            u = u_min;
        double diff = u - u_before_saturation;
        double aw_gain = diff * k_aw;
        prev_aw_gain = aw_gain;

        return { u, e_1 };
    }
}

namespace DT
{
    /*
    PID SYSTEM IMPLEMENTATIONS
    */  
    ClosedLoopSystem_PID::ClosedLoopSystem_PID(DT::TransferFunction *tf, DT::AproximationType aprox_type,
                                               double P, double I, double D, double T, double N,
                                               double uMin, double uMax, double Kaw)
    {
        pid_regulator = std::make_unique<DT::PIDRegulator>(aprox_type, P, I, D, T, N, uMin, uMax, Kaw);
        system = tf;
    }

    ClosedLoopSystem_PID::~ClosedLoopSystem_PID()
    {
        system = nullptr;
    }

    DT::ClosedLoopStepResponse ClosedLoopSystem_PID::step(double w)
    {
        DT::RegulatorResponse reg_res = pid_regulator->step(w, prev_y);
        double y = system->step(reg_res.u);
        prev_y = y;
        return {reg_res.e, reg_res.u, y};
    }

    /*
    PIV SYSTEM IMPLEMENTATIONS
    */
    ClosedLoopSystem_PIV::ClosedLoopSystem_PIV(DT::TransferFunction *tf, DT::AproximationType aprox_type,
                                               double P, double I, double V, double T,
                                               double uMin, double uMax, double Kaw)
    {
        output_integrator = std::make_unique<Integrator>(aprox_type, T);
        piv_regulator = std::make_unique<DT::PIVRegulator>(aprox_type, P, I, V, T, uMin, uMax, Kaw);
        system = tf;
    }

    ClosedLoopSystem_PIV::~ClosedLoopSystem_PIV()
    {
        system = nullptr;
    }

    DT::ClosedLoopStepResponse ClosedLoopSystem_PIV::step(double w)
    {
        DT::RegulatorResponse reg_res = piv_regulator->step(w, prev_y, prev_iy);
        double y = system->step(reg_res.u);             // y  - speed
        double iy = output_integrator->step(y);         // iy - position
        prev_y = y;
        prev_iy = iy;
        return {reg_res.e, reg_res.u, iy};
    }
}