#pragma once

#include "./transfer_fcn.h"
#include <memory>

/*---------- Example of usage PID system regulation in main.cpp: ------------

const Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");

int main(int argc, char const *argv[])
{
    std::ifstream is("input_signal.txt");
    std::istream_iterator<double> start(is), end;
    std::vector<double> in(start, end);
    std::cout << "Read " << in.size() << " numbers" << std::endl;

    Eigen::Vector3d A {{ 1.0, -1.8955, 0.9050 }};
    Eigen::Vector3d B {{ 0.0024, 0.0048, 0.0024 }};
    DT::TransferFunction tf(B, A);

    DT::ClosedLoopSystem_PID cls(&reg, DT::TPZ, 0.7, 0.8, 0.7, 0.1, 10);

    for (int i=0; i<500; i++)
    {
        double w = in[i];
        DT::ClosedLoopStepResponse out = cls.step(w);
        std::cout << out.y << std::endl;
    }

    return 0;
}

/*---------- Example of usage PIV system regulation in main.cpp: ------------

 const Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");

int main(int argc, char const *argv[])
{
    double T_step = 0.1;

    // input
    Eigen::VectorXd in = Eigen::VectorXd::Ones(500);

    // DC motor discrete fcn
    Eigen::VectorXd A {{ 1.0, -0.9802 }};
    Eigen::VectorXd B {{ 0.0594 }};
    DT::TransferFunction discrete_dc_model(B, A);

    Eigen::VectorXd cA {{ 5, 1 }};
    Eigen::VectorXd cB {{ 3 }};
    DT::TransferFunction continuous_dc_model(cB, cA);

    // make pole-placement
    auto PIV = DT::PolePlacement::PIV(continuous_dc_model, 2.0, 0.7, 1.0);  // omega = 2.0, b = 0.7, k = 1.0

    // make PIV closed loop system
    DT::ClosedLoopSystem_PIV cls_piv(&discrete_dc_model, DT::TPZ, PIV.P, PIV.I, PIV.V, T_step);

    std::cout << "Found PIV params: P: " << PIV.P << ", I: " << PIV.I << ", V: " << PIV.V << std::endl;

    // P+IV regulator
    for (int i=0; i<499; i++)
    {
        double w = in[i];
        // double y = discrete_dc_model.step(w);
        DT::ClosedLoopStepResponse res = cls_piv.step(w);
        std::cout << res.y << std::endl; 
    }

    return 0;
}
 
*/

// structures and custom data types
namespace DT
{
    struct ClosedLoopStepResponse
    {
        double e;
        double u;
        double y;
    };

    enum AproximationType
    {
        FWD = 0,
        BWD = 1,
        TPZ = 2,
        PSD = 3
    };
}

// functionality declaration
namespace DT
{

    class Integrator
    {
    private:
        std::unique_ptr<TransferFunction> tf;
    public:
        Integrator(AproximationType type, double T=0.0);
        ~Integrator();

        inline double produceOutput(double in) { return tf->step(in); };
    };

    class Derivator
    {
    private:
        std::unique_ptr<TransferFunction> tf;
    public:
        Derivator(AproximationType type, double T=0.0, double N=0.0);
        ~Derivator();

        inline double produceOutput(double in) { return tf->step(in); };
    };

    class PIDRegulator 
    {
    private:
        double P_gain;
        double I_gain;
        double D_gain;
        std::unique_ptr<Integrator> integrator;
        std::unique_ptr<Derivator> derivator;

    public:
        PIDRegulator(AproximationType aproxType, double P, double I, double D, double T=0.0, double N=0.0);
        ~PIDRegulator();

        double produceOutput(double e);       
    };

    class ClosedLoopSystem_PID
    {
    private:
        std::unique_ptr<DT::PIDRegulator> pid_regulator;
        DT::TransferFunction* system; 
        double previous_y = 0;

    public:
        ClosedLoopSystem_PID(DT::TransferFunction* tf, AproximationType aproxType, double P, double I, double D, double T=0.0, double N=0.0);
        ~ClosedLoopSystem_PID();
        ClosedLoopStepResponse step(double w);
    };

    class ClosedLoopSystem_PIV
    {
    private:
        double P_gain, I_gain, V_gain;
        std::unique_ptr<DT::Integrator> position_integrator;
        std::unique_ptr<DT::Integrator> i_reg_integrator;
        DT::TransferFunction* system; 
        double previous_y = 0;
        double previous_iy = 0;

    public:
        ClosedLoopSystem_PIV(DT::TransferFunction* tf, AproximationType aproxType, double P, double I, double V, double T=0.0);
        ~ClosedLoopSystem_PIV();
        ClosedLoopStepResponse step(double w);
    };
}
