#pragma once

#include "./Eigen/Dense"

namespace DT
{
    class CircleBuffer
    {
    private:
        Eigen::VectorXd vector;
        uint idx;

    public:
        CircleBuffer(Eigen::VectorXd data)
        {
            vector = data;
            idx = vector.size()-1;
        };

        void add(double value) 
        {
            idx++;
            idx%=vector.size();
            vector[idx] = value;           
        };

        double at(uint index)
        {
            if (index > (uint)vector.size()-1) 
            {
#ifdef __x86_64__
                throw std::out_of_range ("Index is out of range!");
#else
                return MAXFLOAT;
#endif
            }
            uint real_index =  (idx+(vector.size()-index))%vector.size();
            return vector[real_index];
        };
    };
}
