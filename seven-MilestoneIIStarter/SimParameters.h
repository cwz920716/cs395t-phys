#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

struct SimParameters
{
    SimParameters()
    {
        timeStep = 0.001;
        NewtonMaxIters = 20;
        NewtonTolerance = 1e-8;
        
        gravityEnabled = true;
        gravityG = 1.0;                
        penaltyEnabled = true;
        penaltyStiffness = 1000.0;
        impulsesEnabled = true;
        CoR = 0.5;

        destroyIndex = 0;
        destroyPieces = 2;
        explosionMag = 10.0;
    }

    float timeStep;
    float NewtonTolerance;
    int NewtonMaxIters;
    
    bool gravityEnabled;
    float gravityG;    
    bool penaltyEnabled;
    float penaltyStiffness;
    bool impulsesEnabled;
    float CoR;

    int destroyIndex;
    int destroyPieces;
    float explosionMag;
};

#endif