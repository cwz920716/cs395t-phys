#ifndef SCENEOBJECTS_H
#define SCENEOBJECTS_H

#include <Eigen/Core>
#include "SimParameters.h"
#include <set>

struct Particle
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Particle(Eigen::Vector2d pos, double mass, bool isFixed, bool isInert) : pos(pos), mass(mass), fixed(isFixed), inert(isInert)
    {
        vel.setZero();
        prevpos = pos;
    }

    Eigen::Vector2d pos;
    Eigen::Vector2d prevpos;
    Eigen::Vector2d vel;
    double mass;
    bool fixed;
    bool inert;
};

struct Connector
{
public:
    Connector(int p1,int p2, double mass) : p1(p1), p2(p2), mass(mass) {}
    virtual ~Connector() {}

    int p1;
    int p2;
    double mass;

    std::set<int> associatedBendingStencils;
};

struct Spring : public Connector
{
public:
    Spring(int p1, int p2, double mass, double stiffness, double restlen, bool canSnap) : Connector(p1, p2, mass), stiffness(stiffness), restlen(restlen), canSnap(canSnap) {}

    double stiffness;
    double restlen;
    bool canSnap;
};

struct Saw
{
public:
    Saw(Eigen::Vector2d pos, double radius) : pos(pos), radius(radius) {}

    Eigen::Vector2d pos;
    double radius;
};

#endif