#ifndef RIGIDBODYTEMPLATE_H
#define RIGIDBODYTEMPLATE_H

#include <string>
#include <Eigen/Core>
#include <set>

class SignedDistanceField;

class RigidBodyTemplate
{
public:
    RigidBodyTemplate(const std::string &meshFilename, double scale);
    ~RigidBodyTemplate();

    double getVolume() const {return volume_;}
    const Eigen::Matrix3d getInertiaTensor() const {return inertiaTensor_;}    

    double getBoundingRadius() const {return radius_;}
    const Eigen::MatrixX3d &getVerts() const {return V;}
    const Eigen::MatrixX3i &getFaces() const {return F;}       

private:
    RigidBodyTemplate(const RigidBodyTemplate &other) = delete;
    RigidBodyTemplate &operator=(const RigidBodyTemplate &other) = delete;

    void computeNorms();
    void computeCOM();
    void computeIntegrals();
    void computeFaceIntegrals(int f);
    void computeProjectionIntegrals(int f);
    void computeCOMProjectionIntegrals(int f);

    void initialize();

    Eigen::MatrixX3d V;
    Eigen::MatrixX3i F;
    
    double volume_;
    double radius_;
    Eigen::Matrix3d inertiaTensor_;

    Eigen::Vector3d c;
    Eigen::MatrixX3d norms;
    Eigen::VectorXd w;
    int numVerts, numFaces;

    int A;   /* alpha */
    int B;   /* beta */
    int C;   /* gamma */

    /* projection integrals */
    double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;

    /* face integrals */
    double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;

    /* volume integrals */
    double T0, T1[3], T2[3], TP[3];
};

#endif // RIGIDBODYTEMPLATE_H
