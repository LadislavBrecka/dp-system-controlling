#include <iostream>

#include "../inc/identification.h"
// #include "../inc/eigen_formatter.h"    /we will use custom one with bigger precious

const Eigen::IOFormat fmt(17, 0, ", ", "\n", "[", "]");

namespace DT 
{

    namespace Examples
    {

        void example_identification() 
        {
            DT::Identificator identificator(DT::LSM, 2, 2);

            Eigen::VectorXd u { { 0.2, 2.5, -3.3, -7.2, 2.8, 5.2, -2.4, 4.0, 8.9, 2.0 } };
            Eigen::VectorXd y { { 0, 0.104600000000000, 1.291737000000000, -1.969669145000000, -4.111832423975001, 
            2.529228559414375, 4.600826351453026, -1.448895839452323, 0.522951885421992, 3.923257076094510 } };

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