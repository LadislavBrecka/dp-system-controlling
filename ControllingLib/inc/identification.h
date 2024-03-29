#pragma once

#include <memory>

#include "./Eigen/Dense"
#include "./common_types.h"

namespace DT 
{

    // main class of identification -> we need instance of this class for identification
    class Identificator;  

    // base class for identification method, which is created in main class and it's this 
    // class which is providing identification (throught it's implementation in child classes)
    class IdentificationMethod;    

    // child classes, which are implementing specific algorithm of identification
    class LeastSquareMethod;       
}

namespace DT 
{
    class IdentificationMethod 
    {
    protected:
        uint n_parameters;
        Eigen::VectorXd thetas;

    public:
        IdentificationMethod(uint n_params);
        ~IdentificationMethod();

        // each child class (each specific method) must implement own update method
        virtual void update(Eigen::VectorXd h, double y) = 0;

        // encapsulation methods
        inline Eigen::VectorXd get_thetas() { return thetas; }
    };

    class Identificator 
    {
    private:
        uint n_a, n_b;
        Eigen::VectorXd h;

        // identification is done by identification method
        // we can choose, which method we want to use
        std::unique_ptr<DT::IdentificationMethod> method;

    public:
        Identificator(DT::IdentificationMethodType method_type, uint num_order, uint den_order);
        ~Identificator();

        // method for updating coeficients of identification
        // we call this funcion in loop and in every call we provide the newest IO data from system
        void update_coeficients(double u, double y);

        // encapsulation methods
        inline Eigen::VectorXd get_thetas() { return method->get_thetas(); }

    private:
        // shifting vector h with new IO data
        void shift_vector_h(double u, double y);
    };
}


namespace DT {

    class LeastSquareMethod : public IdentificationMethod 
    {

    private:
        Eigen::MatrixXd P;
        Eigen::VectorXd d;
        double e, ro, Q;

    public:
        LeastSquareMethod(uint n_params);
        ~LeastSquareMethod();

        void update(Eigen::VectorXd h, double y) override;
    };
}

