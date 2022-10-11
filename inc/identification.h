#pragma once
#include "./Eigen/Dense"
#include <vector>
#include <array>
#include <memory>

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
        uint8_t n_parameters;
        Eigen::VectorXd thetas;

    public:
        IdentificationMethod(uint8_t n_params);
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

        uint8_t n_a; uint8_t n_b;

        Eigen::VectorXd h;
        std::vector<Eigen::Vector4d> H;

    public:
        Identificator(IdentificationMethodType type, uint8_t nominator_order, uint8_t denominator_order);
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
        LeastSquareMethod(uint8_t n_params);
        ~LeastSquareMethod();

        void update(Eigen::VectorXd h, double y) override;
    };
}

