#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

struct SimParameters
{
    SimParameters()
    {
        constraintHandling = CH_PENALTY;
        timeStep = 0.001;
        NewtonMaxIters = 20;
        NewtonTolerance = 1e-8;
        penaltyStiffness = 1e5;

        gravityEnabled = true;
        gravityG = -9.8;
        springsEnabled = true;
        springStiffness = 100;
        maxSpringStrain = 0.2;
        dampingEnabled = true;
        dampingStiffness = 1.0;
        floorEnabled = true;
        bendingEnabled = true;

        clickMode = CM_ADDPARTICLE;
        connectorType = CT_SPRING;
        particleMass = 1.0;
        maxSpringDist = 0.25;
        particleFixed = false;

        rodDensity = 2;
        rodStretchingStiffness = 100;
        rodBendingStiffness = 0.05;
        rodSegments = 5;

        sawRadius= 0.1;
    }

    enum ClickMode {CM_ADDPARTICLE, CM_ADDSAW};

    enum ConstraintHandling {CH_PENALTY, CH_STEPPROJECT, CH_LAGRANGEMULT};

    enum ConnectorType {CT_SPRING, CT_RIGIDROD, CT_FLEXROD};

    ConstraintHandling constraintHandling;
    double timeStep;
    double NewtonTolerance;
    int NewtonMaxIters;
    double penaltyStiffness;

    bool gravityEnabled;
    double gravityG;
    bool springsEnabled;
    bool bendingEnabled;
    double springStiffness;
    double maxSpringStrain;
    bool floorEnabled;
    bool dampingEnabled;
    double dampingStiffness;

    ClickMode clickMode;
    ConnectorType connectorType;
    double particleMass;
    double maxSpringDist;
    bool particleFixed;
    double sawRadius;

    double rodDensity;
    double rodBendingStiffness;
    double rodStretchingStiffness;
    int rodSegments;
};

#endif