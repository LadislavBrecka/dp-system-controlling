#include <cmath>
#include <string>
#include <iostream>

#include "../inc/identification.h"

namespace DT 
{

    //  IDENTIFICATION METHOD BASE CLASS IMPLEMENTATIONS
    IdentificationMethod::IdentificationMethod(uint n_params)
    {
        n_parameters = n_params;
        thetas = Eigen::VectorXd::Zero(n_parameters);
    }

    IdentificationMethod::~IdentificationMethod()
    {
    }

    //  MAIN IDENTIFICATION CLASS IMPLEMENTATIONS
    Identificator::Identificator(IdentificationMethodType method_type, uint num_order, uint den_order)
    : n_a(den_order), n_b(num_order)
    {
        std::cout << "Initializing identification main class!" << std::endl;
        switch (method_type) 
        {
            case LSM:
                method = std::make_unique<LeastSquareMethod>(n_a + n_b);
                break;
            default:
                method = std::make_unique<LeastSquareMethod>(n_a + n_b);
                break;
        }

        h = Eigen::VectorXd::Zero(n_a + n_b);
    }

    Identificator::~Identificator()
    {
    }

    void Identificator::update_coeficients(double u, double y)
    {
        this->method->update(h, y);
        shift_vector_h(u, y);
    }

    void Identificator::shift_vector_h(double u, double y)
    {
        // save actual vector h to H matrix
        H.push_back(h);

        // update vector h with inputs of system
        for (size_t i = n_b - 1; i >= 1; i--) h[i] = h[i-1]; 
        h[0] = u;

        // update vector h with outputs of system
        for (size_t i = n_a + n_b - 1; i >= n_b + 1; i--)  h[i] = h[i-1]; 
        h[n_b] = -y;
    }

    //   IDENTIFICATION METHODS CLASSES
    LeastSquareMethod::LeastSquareMethod(uint n_params) : IdentificationMethod(n_params)
    {
        // P matrix initialization with 10^10
        P = Eigen::MatrixXd::Zero(n_parameters, n_parameters);
        for (size_t i = 0; i < n_parameters; i++)
            P(i,i) = 100000000000;

        // d vector initialization
        d = Eigen::VectorXd::Zero(n_parameters);

        // other scalars initialization
        e = 0; ro = 0; Q = 0;

        std::cout << "............................" << std::endl << std::endl;
        std::cout << "Initialization of LSM!" << std::endl;
        std::cout << "Matrix P:\n" << P.format(DT::Formatter::fmt) << std::endl;
        std::cout << "Vector dT: " << d.transpose().format(DT::Formatter::fmt) << std::endl;
        std::cout << "Scalar e: " << e << std::endl;
        std::cout << "Scalar ro: " << ro << std::endl;
        std::cout << "Scalar Q: " << Q << std::endl;
        std::cout << "............................" << std::endl << std::endl;
    }

    LeastSquareMethod::~LeastSquareMethod()
    {  
    }

    void LeastSquareMethod::update(Eigen::VectorXd h, double y)
    {
        auto hT = h.transpose();

        std::cout << "Vector h: " << hT.format(DT::Formatter::fmt) << std::endl;
        std::cout << "Output y: " << y << std::endl;

        e = y - hT * thetas;
        d = P * h;
        ro = 1 / (1 + hT * d);

        thetas = thetas + ro * e * d;

        P = P - ro * d * hT * P;

        Q = Q + ro * (e*e);
        
        std::cout << "New matrix P:\n" << P.format(DT::Formatter::fmt) << std::endl;
        std::cout << "New thetas:\n" << thetas.transpose().format(DT::Formatter::fmt) << std::endl;
        std:: cout << "Q: " << Q << std::endl;
        std:: cout << "e: " << e << std::endl;
    }
}