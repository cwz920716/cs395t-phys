#include "RigidBodyTemplate.h"
#include <iostream>
#include <igl/readOBJ.h>
#include <Eigen/Dense>
#include <fstream>
#include <map>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

#define X 0
#define Y 1
#define Z 2

#define SQUARE(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

RigidBodyTemplate::RigidBodyTemplate(const std::string &meshFilename, double scale) : volume_(0), radius_(0)
{
    inertiaTensor_.setZero();

    igl::readOBJ(meshFilename, V, F);

    V *= scale;

    initialize();
}

RigidBodyTemplate::~RigidBodyTemplate()
{    
}

void RigidBodyTemplate::computeCOMProjectionIntegrals(int f)
{
    P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0;

    for (int i = 0; i < 3; i++) {
        int ii = (i + 1) % 3;
        double a0 = V(F(f, i), A);
        double b0 = V(F(f, i), B);
        double a1 = V(F(f, ii), A);
        double b1 = V(F(f, ii), B);
        double da = a1 - a0;
        double db = b1 - b0;
        double a0_2 = a0 * a0;
        double a0_3 = a0_2 * a0;
        double b0_2 = b0 * b0;
        double b0_3 = b0_2 * b0;
        double a1_2 = a1 * a1; 
        double b1_2 = b1 * b1;

        double C1 = a1 + a0;
        double Ca = a1*C1 + a0_2;
        double Caa = a1*Ca + a0_3;
        double Cb = b1*(b1 + b0) + b0_2;
        double Cbb = b1*Cb + b0_3;
        double Cab = 3*a1_2 + 2*a1*a0 + a0_2;
        double Kab = a1_2 + 2*a1*a0 + 3*a0_2;

        P1 += db*C1;
        Pa += db*Ca;
        Paa += db*Caa;
        Pb += da*Cb;
        Pbb += da*Cbb;
        Pab += db*(b1*Cab + b0*Kab);
    }

    P1 /= 2.0;
    Pa /= 6.0;
    Paa /= 12.0;
    Paaa /= 20.0;
    Pb /= -6.0;
    Pbb /= -12.0;
    Pbbb /= -20.0;
    Pab /= 24.0;
    Paab /= 60.0;
    Pabb /= -60.0;
}

void RigidBodyTemplate::computeProjectionIntegrals(int f)
{
    P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0;

    for (int i = 0; i < 3; i++) {
        int ii = (i + 1) % 3;
        double a0 = V(F(f, i), A);
        double b0 = V(F(f, i), B);
        double a1 = V(F(f, ii), A);
        double b1 = V(F(f, ii), B);
        double da = a1 - a0;
        double db = b1 - b0;
        double a0_2 = a0 * a0;
        double a0_3 = a0_2 * a0;
        double a0_4 = a0_3 * a0;
        double b0_2 = b0 * b0;
        double b0_3 = b0_2 * b0;
        double b0_4 = b0_3 * b0;
        double a1_2 = a1 * a1;
        double a1_3 = a1_2 * a1; 
        double b1_2 = b1 * b1;
        double b1_3 = b1_2 * b1;

        double C1 = a1 + a0;
        double Ca = a1*C1 + a0_2;
        double Caa = a1*Ca + a0_3;
        double Caaa = a1*Caa + a0_4;
        double Cb = b1*(b1 + b0) + b0_2;
        double Cbb = b1*Cb + b0_3;
        double Cbbb = b1*Cbb + b0_4;
        double Cab = 3*a1_2 + 2*a1*a0 + a0_2;
        double Kab = a1_2 + 2*a1*a0 + 3*a0_2;
        double Caab = a0*Cab + 4*a1_3;
        double Kaab = a1*Kab + 4*a0_3;
        double Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
        double Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

        P1 += db*C1;
        Pa += db*Ca;
        Paa += db*Caa;
        Paaa += db*Caaa;
        Pb += da*Cb;
        Pbb += da*Cbb;
        Pbbb += da*Cbbb;
        Pab += db*(b1*Cab + b0*Kab);
        Paab += db*(b1*Caab + b0*Kaab);
        Pabb += da*(a1*Cabb + a0*Kabb);
    }

    P1 /= 2.0;
    Pa /= 6.0;
    Paa /= 12.0;
    Paaa /= 20.0;
    Pb /= -6.0;
    Pbb /= -12.0;
    Pbbb /= -20.0;
    Pab /= 24.0;
    Paab /= 60.0;
    Pabb /= -60.0;
}

void RigidBodyTemplate::computeFaceIntegrals(int f) {
    double k1 = 0, k2 = 0, k3 = 0, k4 = 0;

    computeProjectionIntegrals(f);
    k1 = 1 / norms(f, C); k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;

    Fa = k1 * Pa;
    Fb = k1 * Pb;
    Fc = -k2 * (norms(f, A)*Pa + norms(f, B)*Pb + w(f)*P1);

    Faa = k1 * Paa;
    Fbb = k1 * Pbb;
    Fcc = k3 * (SQUARE(norms(f, A)) * Paa +
                2 * norms(f, A) * norms(f, B) * Pab +
                SQUARE(norms(f, B)) * Pbb +
                2 * w(f) * norms(f, A) * Pa +
                2 * w(f) * norms(f, B) * Pb +
                SQUARE(w(f)) * P1);

    Faaa = k1 * Paaa;
    Fbbb = k1 * Pbbb;
    Fccc = -k4 * (CUBE(norms(f, A)) * Paaa +
                  3 * SQUARE(norms(f, A)) * norms(f, B) * Paab +
                  3 * norms(f, A) * SQUARE(norms(f, B)) * Pabb +
                  CUBE(norms(f, B)) * Pbbb +
                  3 * SQUARE(norms(f, A)) * w(f) * Paa +
                  3 * SQUARE(norms(f, B)) * w(f) * Pbb +
                  6 * norms(f, A) * w(f) * norms(f, B) * Pab +
                  3 * norms(f, A) * SQUARE(w(f)) * Pa +
                  3 * norms(f, B) * SQUARE(w(f)) * Pb +
                  CUBE(w(f)) * P1);

    Faab = k1 * Paab;
    Fbbc = -k2 * (norms(f, A) * Pabb + norms(f, B) * Pbbb + w(f) * Pbb);
    Fcca = k3 * (SQUARE(norms(f, A)) * Paaa +
                 2 * norms(f, A) * norms(f, B) * Paab +
                 SQUARE(norms(f, B))*Pabb +
                 2 * norms(f, A) * w(f) * Paa +
                 2 * norms(f, B) * w(f) * Pab +
                 SQUARE(w(f)) * Pa);
}

void RigidBodyTemplate::computeNorms() {
    for (int i = 0; i < numFaces; i++) {

        // std::cout << "v0 = [" << V(F(i, 0), X) << "," << V(F(i, 0), Y) << "," << V(F(i, 0), Z) << "]\n";
        // std::cout << "v1 = [" << V(F(i, 1), X) << "," << V(F(i, 1), Y) << "," << V(F(i, 1), Z) << "]\n";
        // std::cout << "v2 = [" << V(F(i, 2), X) << "," << V(F(i, 2), Y) << "," << V(F(i, 2), Z) << "]\n";

        double dx1 = V(F(i, 1), X) - V(F(i, 0), X);
        double dy1 = V(F(i, 1), Y) - V(F(i, 0), Y);
        double dz1 = V(F(i, 1), Z) - V(F(i, 0), Z);
        double dx2 = V(F(i, 2), X) - V(F(i, 1), X);
        double dy2 = V(F(i, 2), Y) - V(F(i, 1), Y);
        double dz2 = V(F(i, 2), Z) - V(F(i, 1), Z);
        double nx = dy1 * dz2 - dy2 * dz1;
        double ny = dz1 * dx2 - dz2 * dx1;
        double nz = dx1 * dy2 - dx2 * dy1;
        double len = sqrt(nx * nx + ny * ny + nz * nz);
        norms(i, X) = nx / len;
        norms(i, Y) = ny / len;
        norms(i, Z) = nz / len;

        w(i) = -(norms(i, X) * V(F(i, 0), X) + 
                 norms(i, Y) * V(F(i, 0), Y) +
                 norms(i, Z) * V(F(i, 0), Z));

        // std::cout << "norm = " << norms.block(i, 0, 1, 3) << "\nw = " << w(i) << "\n";
    }
}

void RigidBodyTemplate::computeCOM() {
    T0 = T1[X] = T1[Y] = T1[Z] = 0;

    for (int i = 0; i < numFaces; i++) {

        double nx = fabs(norms(i, X));
        double ny = fabs(norms(i, Y));
        double nz = fabs(norms(i, Z));
        if (nx > ny && nx > nz) C = X;
        else C = (ny > nz) ? Y : Z;
        A = (C + 1) % 3;
        B = (A + 1) % 3;

        double k1 = 0, k2 = 0, k3 = 0;

        computeCOMProjectionIntegrals(i);
        k1 = 1 / norms(i, C); k2 = k1 * k1; k3 = k2 * k1;

        Fa = k1 * Pa;
        Fb = k1 * Pb;
        Fc = -k2 * (norms(i, A) * Pa + norms(i, B) * Pb + w(i) * P1);

        Faa = k1 * Paa;
        Fbb = k1 * Pbb;
        Fcc = k3 * (SQUARE(norms(i, A)) * Paa +
                    2 * norms(i, A) * norms(i, B) * Pab +
                    SQUARE(norms(i, B)) * Pbb +
                    2 * w(i) * norms(i, A) * Pa +
                    2 * w(i) * norms(i, B) * Pb +
                    SQUARE(w(i)) * P1);

        if (A == X) {
            T0 += norms(i, A) * Fa;
        } else if (B == X) {
            T0 += norms(i, B) * Fb;
        } else {
            T0 += norms(i, C) * Fc;
        }

        T1[A] += norms(i, A) * Faa;
        T1[B] += norms(i, B) * Fbb;
        T1[C] += norms(i, C) * Fcc;
    }

    T1[X] /= 2; T1[Y] /= 2; T1[Z] /= 2;
}

void RigidBodyTemplate::computeIntegrals() {

    T0 = T1[X] = T1[Y] = T1[Z] 
       = T2[X] = T2[Y] = T2[Z] 
       = TP[X] = TP[Y] = TP[Z] = 0;

    for (int i = 0; i < numFaces; i++) {

        double nx = fabs(norms(i, X));
        double ny = fabs(norms(i, Y));
        double nz = fabs(norms(i, Z));
        if (nx > ny && nx > nz) C = X;
        else C = (ny > nz) ? Y : Z;
        A = (C + 1) % 3;
        B = (A + 1) % 3;

        computeFaceIntegrals(i);

        if (A == X) {
            T0 += norms(i, A) * Fa;
        } else if (B == X) {
            T0 += norms(i, B) * Fb;
        } else {
            T0 += norms(i, C) * Fc;
        }

        T1[A] += norms(i, A) * Faa;
        T1[B] += norms(i, B) * Fbb;
        T1[C] += norms(i, C) * Fcc;
        // Taa
        T2[A] += norms(i, A) * Faaa;
        // Tbb
        T2[B] += norms(i, B) * Fbbb;
        // Tcc
        T2[C] += norms(i, C) * Fccc;
        // Tab
        TP[A] += norms(i, A) * Faab;
        // Tbc
        TP[B] += norms(i, B) * Fbbc;
        // Tac
        TP[C] += norms(i, C) * Fcca;
    }

    T1[X] /= 2; T1[Y] /= 2; T1[Z] /= 2;
    T2[X] /= 3; T2[Y] /= 3; T2[Z] /= 3;
    TP[X] /= 2; TP[Y] /= 2; TP[Z] /= 2;
}

void RigidBodyTemplate::initialize()
{
  // TODO compute quantities such as volume, center of mass, and intertia tensor
  // translate body so its center of mass is the world origin

    numVerts = V.rows();
    numFaces = F.rows();
    Eigen::MatrixX3d m(numFaces, 3);
    Eigen::VectorXd v(numFaces);
    norms = m;
    w = v;

    // std::cout << "V=" << numVerts << ", F=" << numFaces << "\n";
    computeNorms();
    computeCOM();

    volume_ = T0;
    // std::cout << "V=" << volume_ << "\n";

    c(X) = T1[X] / volume_;
    c(Y) = T1[Y] / volume_;
    c(Z) = T1[Z] / volume_;
    // std::cout << "C=[\n" << c << "]\n";

    for (int i = 0; i < numVerts; i++) {
        V(i, X) -= c(X);
        V(i, Y) -= c(Y);
        V(i, Z) -= c(Z);
    }
    computeNorms();
    computeIntegrals();

    /* compute inertia tensor */
    inertiaTensor_(X, X) = T2[Y] + T2[Z];
    inertiaTensor_(Y, Y) = T2[Z] + T2[X];
    inertiaTensor_(Z, Z) = T2[X] + T2[Y];
    inertiaTensor_(X, Y) = inertiaTensor_(Y, X) = - TP[X];
    inertiaTensor_(Y, Z) = inertiaTensor_(Z, Y) = - TP[Y];
    inertiaTensor_(Z, X) = inertiaTensor_(X, Z) = - TP[Z];

    // std::cout << "inertiaTensor_=\n" << inertiaTensor_ << "\n";
}    
