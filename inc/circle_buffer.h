#pragma once
#include "Eigen/Dense"

namespace DT{
    class CircleBuffer {

    private:
        Eigen::VectorXd vector;
        uint size;
        uint idx;

    public:
        CircleBuffer(Eigen::VectorXd data)
        {
            vector = data;
            size = data.size();
            idx = UINT_MAX;
        };

        void add(double value) 
        {
            if (idx >= size-1) idx = 0; else idx++;                      
            vector[idx] = value;           
        };

        double at(uint index)
        {
            if (index > size-1) throw std::out_of_range ("Index is out of range!");
            int real_index = idx - index;
            if (real_index < 0) real_index = size + real_index;
            return vector[real_index];
        }

        // double &operator[] (uint index) {
        //     if (index > size-1) throw std::out_of_range ("Index is out of range!");

        //     int real_index = idx - index;
        //     if (real_index < 0) real_index = size + real_index;
        //     return vector[real_index];
        // };
    };
}
