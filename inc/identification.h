#pragma once
#include "./Eigen/Dense"
#include <vector>
#include <array>
#include <memory>

/* 
---------- Example of usage in main.cpp: ------------

const Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");

int main(int argc, char const *argv[])
{
    DT::Identificator identificator(DT::LSM, 2, 2);

    Eigen::VectorXd u { {  1.0,  0.0,   -5.5,     1.5,      3.3,     10.0,   -1.0,     1.1,     2.5,   -4.8   } };
    Eigen::VectorXd y { {   0,  0.006,  0.0154, -0.0093,  -0.0445,  -0.0502, 0.0184,  0.1143,  0.2025, 0.2998 } };

    for (int i = 0; i < u.size(); i++)
    {
        std::cout << "\n-------- New iteration " << i+1 << " ---------" << std::endl;
        identificator.updateCoeficients(u[i], y[i]);
    }

    auto thetas = identificator.getThetas();

    std::cout << std::endl << "----------------------------------------" << std::endl;
    std::cout << std::endl << "|          Found parameters:           |" << std::endl;
    std::cout << std::endl << "----------------------------------------" << std::endl;
    std::cout << thetas.format(fmt);

    return 0;
}

Compile command: g++ main.cpp src/identification.cpp

*/


namespace DT {

    // main class of identification -> we need instance of this class for identification
    class Identificator;  

    // base class for identification method, which is created in main class and it's this 
    // class which is providing identification (throught it's implementation in child classes)
    class IdentificationMethod;    

    // child classes, which are implementing specific algorithm of identification
    class LeastSquareMethod;       

    // enum just for method types distinguishing 
    enum IdentificationMethodType { LSM = 0 };

    // ----------------------------------- //
    //  IDENTIFICATION METHOD BASE CLASS  //
    // ----------------------------------- //
    class IdentificationMethod 
    {

    protected:
        uint n_parameters;
        Eigen::VectorXd thetas;

    public:
        IdentificationMethod(uint n_params);
        ~IdentificationMethod();

        // each child class (each specific method must implement own update methods)
        virtual void update(Eigen::VectorXd h, double y) = 0;

        // encapsulation methods
        inline Eigen::VectorXd getThetas() { return thetas; }
    };

    // --------------------------- //
    //  MAIN IDENTIFICATION CLASS  //
    // --------------------------- //
    class Identificator 
    {

    private:
        // identification is done by identification method
        // we can choose, which method we want to use
        std::unique_ptr<IdentificationMethod> method;

        uint n_a; uint n_b;

        Eigen::VectorXd h;
        std::vector<Eigen::Vector4d> H;

    public:
        Identificator(IdentificationMethodType type, uint nominator_order, uint denominator_order);
        ~Identificator();

        // method for updating coeficients of identification
        // we call this funcion in loop and in every call we provide the newest IO data from system
        void updateCoeficients(double u, double y);

        // encapsulation methods
        inline Eigen::VectorXd getThetas() { return method->getThetas(); }
        inline std::vector<Eigen::Vector4d> getMatrixH() { return H; }

    private:
        // shifting vector h with new IO data
        void shiftVector_h(double u, double y);
    };


    // ----------------------------------- //
    //  IDENTIFICATION METHODS CLASSES     //
    // ----------------------------------- // 

    // LEAST SQUARE METHOD
    class LeastSquareMethod : public IdentificationMethod 
    {

    private:
        Eigen::MatrixXd P;
        Eigen::VectorXd d;
        double e;
        double ro;
        double Q;

    public:
        LeastSquareMethod(uint n_params);
        ~LeastSquareMethod();

        void update(Eigen::VectorXd h, double y) override;
    };
}

