#pragma once
#include "Eigen/Dense"

namespace DT{
    class CircleBuffer {

    private:
        Eigen::VectorXd vector;
        // uint size;
        uint idx;

    public:
        CircleBuffer(Eigen::VectorXd data)
        {
            vector = data;
            // size = data.size();  skusme to takto - nevyrabat informacie ktore su obsiahnute v inom objekte
            idx = vector.size()-1;  //UINT_MAX; nemalo by sa stat nic zle ak to dam takto
        };

        void add(double value) 
        {
            //if (idx >= vector.size()-1) idx = 0; else idx++;
            idx++;
            idx%=vector.size(); //skuste tento trik s modulo operatorom
            vector[idx] = value;           
        };

        double at(uint index)
        {
            if (index > vector.size()-1) throw std::out_of_range ("Index is out of range!");
            int real_index = idx - index;
            if (real_index < 0) real_index = vector.size() + real_index; 
            return vector[real_index];
        }
    };
}
