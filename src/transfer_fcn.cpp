#include <cmath>
#include <iostream>

#include "../inc/transfer_fcn.h"
#include "../inc/Exceptions/not_supported_exception.h"
#include "../inc/Eigen/unsupported/Polynomials"

std::string convert_to_uper(int index, std::string var);

namespace DT {

    TransferFunction::TransferFunction()
    {
        set_numerator(Eigen::VectorXd::Ones(1));
        set_denominator(Eigen::VectorXd::Ones(1));
    }

    TransferFunction::TransferFunction(Eigen::VectorXd num, Eigen::VectorXd den)
    {
        set_numerator(num);
        set_denominator(den);
    }
    
    TransferFunction::TransferFunction(std::vector<double> num, std::vector<double> den)
    {
        set_numerator(num);
        set_denominator(den);
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
            y = (input_part - output_part)/A[0];
          
        // shift vector y so it containts the newest output sample
        vY->add(y);
        return y;
    }
    
    // TODO: not properly tested!!!
    void TransferFunction::d2c(double Ts, DT::TransferFunction& c_tf)
    {
        // fill companion matrix
        Eigen::MatrixXd companion_matrix = Eigen::MatrixXd::Zero(n_a-1, n_a-1);  
        for (uint i=0; i<n_a-1; i++)        
        {
             companion_matrix(i, n_a - 2) = -A[n_a - 1 - i];  
             if(i+1<n_a-1) companion_matrix(i+1, i) = 1;
        }
        // find eigenvalue of companion matrix -> found values are roots of polynom
        Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(companion_matrix);
        Eigen::VectorXcd roots = eigensolver.eigenvalues();
        
        // convert discrete roots to continuous roots
        Eigen::VectorXcd c_roots = Eigen::VectorXcd::Zero(roots.size());
        for (uint i = n_a - 1; i > 0; i--)
        {
            std::complex<double> log_value(log(std::abs(roots(i-1))), std::arg(roots(i-1))); 
            c_roots(i-1) = log_value / Ts;
        }

        // make polynom from continuous roots
        Eigen::VectorXcd c_poly(n_a);
        Eigen::roots_to_monicPolynomial(c_roots, c_poly);

        Eigen::VectorXd c_A = c_poly.real().reverse();
        Eigen::VectorXd c_B {{ B.sum() / A.sum() }};

        // set a_0 as 1
        c_A = c_A /c_A(n_a-1);  

        c_tf.set_denominator(c_A);
        c_tf.set_numerator(c_B);
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
    
    void TransferFunction::set_numerator(Eigen::VectorXd numerator)
    {
        n_b = numerator.size();
        B = numerator;
        vU = std::make_unique<DT::CircleBuffer>(Eigen::VectorXd::Zero(n_b));
    }
    
    void TransferFunction::set_numerator(std::vector<double> numerator)
    {
        n_b = numerator.size();
        B = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(numerator.data(), n_b);
        vU = std::make_unique<DT::CircleBuffer>(Eigen::VectorXd::Zero(n_b));
    }
    
    void TransferFunction::set_denominator(Eigen::VectorXd denominator)
    {
        n_a = denominator.size();
        A = denominator;
        vY = std::make_unique<DT::CircleBuffer>(Eigen::VectorXd::Zero(n_a));
    }
    
    void TransferFunction::set_denominator(std::vector<double> denominator)
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
