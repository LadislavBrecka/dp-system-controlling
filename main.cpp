#include "inc/identification.h"
#include "inc/transfer_fcn.h"
#include <iostream>

const Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");

int main(int argc, char const *argv[])
{
    DT::Identificator identificator(DT::LSM, 2, 2);

    Eigen::Vector3d A {{ 1, -1.895, 0.9048 }};
    Eigen::Vector3d B {{ 0, 0.006, 0.004 }};
    DT::TransferFunction tf(B, A);

    Eigen::VectorXd u { {  1.0,  0.0,   -5.5,     1.5,      3.3,     10.0,   -1.0,     1.1,     2.5,   -4.8   } };

    for (int i = 0; i < u.size(); i++)
    {
        std::cout << "\n-------- New iteration " << i+1 << " ---------" << std::endl;
        double y = tf.step(u[i]);
        identificator.updateCoeficients(u[i], y);
    }

    auto thetas = identificator.getThetas();

    std::cout << std::endl << "----------------------------------------" << std::endl;
    std::cout << std::endl << "|          Found parameters:           |" << std::endl;
    std::cout << std::endl << "----------------------------------------" << std::endl;
    std::cout << thetas.format(fmt) << std::endl;

    return 0;
}
