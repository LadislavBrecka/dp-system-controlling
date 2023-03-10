#include <iostream>

#include "../inc/identification.h"
#include "../inc/eigen_formatter.h"

namespace DT 
{

    namespace Examples
    {

        void example_identification() 
        {
            DT::Identificator identificator(DT::LSM, 2, 2);

            Eigen::VectorXd u { {  1.0,  0.0,   -5.5,     1.5,      3.3,     10.0,   -1.0,     1.1,     2.5,   -4.8   } };
            Eigen::VectorXd y { {   0,  0.006,  0.0154, -0.0093,  -0.0445,  -0.0502, 0.0184,  0.1143,  0.2025, 0.2998 } };

            for (int i = 0; i < u.size(); i++)
            {
                std::cout << "\n-------- New iteration " << i+1 << " ---------" << std::endl;
                identificator.update_coeficients(u[i], y[i]);
            }

            auto thetas = identificator.get_thetas();

            std::cout << std::endl << "----------------------------------------" << std::endl;
            std::cout << std::endl << "|          Found parameters:           |" << std::endl;
            std::cout << std::endl << "----------------------------------------" << std::endl;
            std::cout << thetas.format(fmt) << std::endl;
        }
    }
}