#include "../inc/transfer_fcn.h"
#include <iostream>
#include <cmath>

namespace DT {

    TransferFunction::TransferFunction()
    {
        setNominator(Eigen::VectorXd::Ones(1));
        setDenominator(Eigen::VectorXd::Ones(1));
    }

    TransferFunction::TransferFunction(Eigen::VectorXd nominator, Eigen::VectorXd denominator)
    {
        setNominator(nominator);
        setDenominator(denominator);
    }
    
    TransferFunction::TransferFunction(std::vector<double> nominator, std::vector<double> denominator)
    {
        setNominator(nominator);
        setDenominator(denominator);
    }
    
    TransferFunction::~TransferFunction()
    {
    }
    
    double TransferFunction::step(double u)
    {
        // shift vector u so it contains the newest input sample
        vU->add(u);

        // input part
        double input_part = 0;
        for (uint i = 0; i < n_b; i++)
        {
            input_part += B[i] * vU->at(i);
        }

        double output_part = 0;
        for (uint i = 1; i < n_a; i++)
        {
            output_part += A[i] * vY->at(i-1);
        }

        double y = 0.0;
        if (A[0] != 0.0)
            y = (1 / A[0]) * (input_part - output_part);
          
        // shift vector y so it containts the newest output sample
        vY->add(y);
        return y;
    }
    
    // TODO: not properly tested!!!
    void TransferFunction::d2c(double Ts, DT::TransferFunction& c_tf)
    {
        double A_sum = A.sum();
        double B_sum = B.sum();
        double gain = B_sum/A_sum;

        // fill companion matrix
        Eigen::MatrixXd companion_matrix = Eigen::MatrixXd::Zero(n_a-1, n_a-1);
        for (uint i=0; i<n_a-1; i++)   // cols
        {
            for (uint j=0; j<n_a-1; j++)   // rows
            {
                if (i == n_a - 2)
                    companion_matrix(j, i) = -A[n_a - 1 - j];

                if (i+1 == j)
                    companion_matrix(j, i) = 1;
            }
        }
        std::cout << companion_matrix << std::endl;

        // find eigenvalue of companion matrix -> found values are roots of polynom
        Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(companion_matrix);
        Eigen::VectorXcd roots = eigensolver.eigenvalues();
        std::cout << roots << std::endl;

        // Using Eigen unsupported module - example of alternative way of finding roots with Eigen
        // Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
        // solver.compute(A.reverse());
        // const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &r = solver.roots();
        // std::cout << r << std::endl;

        Eigen::VectorXd continuous_roots(n_a - 1);
        Eigen::VectorXd c_A = Eigen::VectorXd::Zero(n_a);
        Eigen::VectorXd c_B {{ gain }};
        double multiplier;

        for (int i = n_a - 1; i > 0; i--)
        {
            double temp = log(roots(i-1).real()) / Ts;
            std::cout << temp << std::endl;
            if (i == n_a - 1)
            {
                multiplier = 1.0 / temp;
                
            }

            c_A(i) = temp * multiplier;
               
            std::cout << i << std::endl;
        }

        c_A(0) = -multiplier;

        std::cout << "Gain of TF: "<< gain << std::endl;
        std::cout << "Roots of TF: " << continuous_roots * multiplier << std::endl;

        c_tf.setDenominator(c_A);
        c_tf.setNominator(c_B);
        c_tf.print();
    }
    
    void TransferFunction::print()
    {
        // nominator printing
        for (uint i = 0; i < n_b; i++)
        {
            if (B[i] != 0.0)
            {
                if (B[i] > 0.0) std::cout << " + "; else std::cout << " - ";
                std::string z_index = (i == 0 ? "" : "z-" + std::to_string(i));
                std::cout << fabs(B[i]) << z_index;
            }
        }

        // dividing line printing
        std::cout << std::endl;
        for (uint i = 0; i < std::max(n_a, n_b) * 8; i++)
        {
            std::cout << "-";
        }
        std::cout << std::endl;

        // denominator printing
        for (uint i = 0; i < n_a; i++)
        {
            if (A[i] != 0.0)
            {
                if (A[i] > 0.0) std::cout << " + "; else std::cout << " - ";
                std::string z_index = (i == 0 ? "" : "z-" + std::to_string(i));
                std::cout << fabs(A[i]) << z_index;
            }
        }

        // new line at the end
        std::cout << std::endl;
    }
    
    void TransferFunction::setNominator(Eigen::VectorXd nominator)
    {
        n_b = nominator.size();
        B = nominator;
        vU = std::make_unique<DT::CircleBuffer>(Eigen::VectorXd::Zero(n_b));
    }
    
    void TransferFunction::setNominator(std::vector<double> nominator)
    {
        n_b = nominator.size();
        B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(nominator.data(), n_b);
        vU = std::make_unique<DT::CircleBuffer>(Eigen::VectorXd::Zero(n_b));
    }
    
    void TransferFunction::setDenominator(Eigen::VectorXd denominator)
    {
        n_a = denominator.size();
        A = denominator;
        vY = std::make_unique<DT::CircleBuffer>(Eigen::VectorXd::Zero(n_a));
    }
    
    void TransferFunction::setDenominator(std::vector<double> denominator)
    {
        n_a = denominator.size();
        A = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(denominator.data(), n_a);
        vY = std::make_unique<DT::CircleBuffer>(Eigen::VectorXd::Zero(n_a));
    }

}