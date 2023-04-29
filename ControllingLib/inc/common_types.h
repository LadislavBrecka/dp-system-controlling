#pragma once

namespace DT
{
    enum IdentificationMethodType
    { 
        LSM = 0
    };

    enum AproximationType
    {
        TPZ = 2,
        PSD = 3
    };

    struct ClosedLoopStepResponse
    {
        double e;
        double u;
        double y;
    };

    struct RegulatorResponse
    {
        double u;
        double e;        
    };

        struct PIDRegCoefs
    {
        double P;
        double I;
        double D;
    };

    struct PIVRegCoefs
    {
        double P;
        double I;
        double V;
    };
}
