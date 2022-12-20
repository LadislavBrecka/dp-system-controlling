#include "../inc/transfer_fcn.h"
#include "../inc/Exceptions/not_supported_exception.h"
#include "../inc/unsupported/Eigen/Polynomials"
#include <iostream>
#include <cmath>

std::string convert_to_uper(int index, std::string var);

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
    // TODO: not working with complex roots!!!
    void TransferFunction::d2c(double Ts, DT::TransferFunction& c_tf)
    {
        double A_sum = A.sum();
        double B_sum = B.sum();
        double gain = B_sum / A_sum;

        // fill companion matrix
        Eigen::MatrixXd companion_matrix = Eigen::MatrixXd::Zero(n_a-1, n_a-1);
        for (uint i=0; i<n_a-1; i++)        // cols
        {
            for (uint j=0; j<n_a-1; j++)    // rows
            {
                if (i == n_a - 2)
                    companion_matrix(j, i) = -A[n_a - 1 - j];

                if (i+1 == j)
                    companion_matrix(j, i) = 1;
            }
        }

        // find eigenvalue of companion matrix -> found values are roots of polynom
        Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(companion_matrix);
        Eigen::VectorXcd roots = eigensolver.eigenvalues();

        // Using Eigen unsupported module - example of alternative way of finding roots with Eigen
        // Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
        // solver.compute(A.reverse());
        // const Eigen::PolynomialSolver<double, Eigen::Dynamic>::RootsType &r = solver.roots();
        // std::cout << r << std::endl;

        Eigen::VectorXcd c_poly(n_a);
        Eigen::VectorXcd c_roots = Eigen::VectorXcd::Zero(roots.size());

        // convert discrete roots to continuous roots
        for (uint i = n_a - 1; i > 0; i--)
        {
            // complex natural logaritmus
            double real_part= roots(i-1).real();
            double imag_part = roots(i-1).imag(); 
            double manginute = log(sqrt(pow(real_part, 2) + pow(imag_part, 2)));
            double phase = std::atan2(imag_part,real_part);
            std::complex<double> log_value(manginute, phase);

            c_roots(i-1) = log_value.real() / Ts;              
        }

        // make polynom from continuous roots
        Eigen::roots_to_monicPolynomial(c_roots, c_poly);

        // std::cout << "Discrete roots are: " << roots << std::endl;
        // std::cout << "Continuos roots are: " << c_roots << std::endl;
        // std::cout << "Continuos polynom is: " << c_poly << std::endl;

        Eigen::VectorXd c_A = c_poly.real().reverse();
        Eigen::VectorXd c_B {{ gain }};

        double multiplier = 1.0 / c_A(n_a-1);
        c_A = c_A * multiplier;

        c_tf.setDenominator(c_A);
        c_tf.setNominator(c_B);
    }
    
    void TransferFunction::print(const std::string& var)
    {
        
        // nominator printing
        std::cout << std::endl << std::endl;
        for (uint i = 0; i < n_b; i++)
        {
            if (B[i] != 0.0)
            {
                if (B[i] > 0.0 && i != 0) 
                    std::cout << " + "; 
                else if (B[i] > 0.0 && i == 0 ) 
                    std::cout << "";
                else 
                    std::cout << "-";   

                std::string indexed_variable = convert_to_uper(n_b-1-i, var);
                if (i == 0 && B[i] == 1)
                    std::cout << indexed_variable;
                else                
                    std::cout << fabs(B[i]) << indexed_variable;
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
                if (A[i] > 0.0 && i != 0) 
                    std::cout << " + "; 
                else if (A[i] > 0.0)
                    std::cout << "";
                else 
                    std::cout << " - ";

                std::string indexed_variable = convert_to_uper(n_a-1-i, var);
                if (i == 0 && A[i] == 1)
                    std::cout << indexed_variable;
                else                
                    std::cout << fabs(A[i]) << indexed_variable;
            }
        }

        // new line at the end
        std::cout << std::endl << std::endl;
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

std::string convert_to_uper(int index, std::string var) 
{
    switch(abs(index))
    {
        case 0: return "";
        case 1: return var;
        case 2: return var + "^2";
        case 3: return var + "^3";
        case 4: return var + "^4";
        case 5: return var + "^5";
        default: return "";
    }
}