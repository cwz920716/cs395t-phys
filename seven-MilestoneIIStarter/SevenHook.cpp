#include "SevenHook.h"
#include "RigidBodyTemplate.h"
#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include "igl/opengl/glfw/imgui/ImGuiHelpers.h"
#include "CollisionDetection.h"
#include <Eigen/Geometry>

using namespace Eigen;

void SevenHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputText("Filename", sceneFile_);
    }
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputFloat("Timestep", &params_.timeStep, 0, 0, 3);
        ImGui::DragFloat("Newton Tolerance", &params_.NewtonTolerance, 0.01, 1e-16, 1e-1, "%.3e", 10);
        ImGui::InputInt("Newton Max Iters", &params_.NewtonMaxIters);
    }
    if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Gravity Enabled", &params_.gravityEnabled);
        ImGui::InputFloat("Gravity G", &params_.gravityG, 0, 0, 3);
        ImGui::Checkbox("Penalty Forces Enabled", &params_.penaltyEnabled);
        ImGui::InputFloat("Penalty Stiffness", &params_.penaltyStiffness, 0, 0, 3);
        ImGui::Checkbox("Impulses Enabled", &params_.impulsesEnabled);
        ImGui::InputFloat("CoR", &params_.CoR, 0, 0, 3);
    }
    if (ImGui::CollapsingHeader("Explosion Controls", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputInt("Object", &params_.destroyIndex);
        ImGui::InputInt("Pieces", &params_.destroyPieces);
        ImGui::InputFloat("Magnitude", &params_.explosionMag, 0, 0, 3);
        if (ImGui::Button("Explode!", ImVec2(-1, 0)))
        {
            destroyMutex_.lock();
            destroyCommands_.push_back(params_.destroyIndex);
            destroyMutex_.unlock();            
        }
    }
}

void SevenHook::updateRenderGeometry()
{
    int totverts = 0;
    int totfaces = 0;
    for (RigidBodyInstance *rbi : bodies_)
    {
        totverts += rbi->getTemplate().getVerts().rows();
        totfaces += rbi->getTemplate().getFaces().rows();
    }
    renderQ.resize(totverts, 3);
    renderF.resize(totfaces, 3);
    int voffset = 0;
    int foffset = 0;
    for (RigidBodyInstance *rbi : bodies_)
    {
        int nverts = rbi->getTemplate().getVerts().rows();
        for (int i = 0; i < nverts; i++)
            renderQ.row(voffset + i) = (rbi->c + VectorMath::rotationMatrix(rbi->theta)*rbi->getTemplate().getVerts().row(i).transpose()).transpose();
        int nfaces = rbi->getTemplate().getFaces().rows();
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                renderF(foffset + i, j) = rbi->getTemplate().getFaces()(i, j) + voffset;
            }
        }
        voffset += nverts;
        foffset += nfaces;
    }
}


void SevenHook::initSimulation()
{
    time_ = 0;    
    loadScene();
    updateRenderGeometry();
}

void SevenHook::tick()
{    
    destroyMutex_.lock();
    for(int body : destroyCommands_)        
    { 
        //TODO destroy object "body"
    }
    destroyCommands_.clear();
    destroyMutex_.unlock();
}

void SevenHook::computeForces(VectorXd &Fc, VectorXd &Ftheta)
{
    Fc.resize(3*bodies_.size());
    Ftheta.resize(3*bodies_.size());
    Fc.setZero();
    Ftheta.setZero();    

    if(params_.gravityEnabled)
    {
        for(int i=0; i<bodies_.size(); i++)
        {
            for(int j=i+1; j<bodies_.size(); j++)
            {
                Vector3d diff = bodies_[j]->c - bodies_[i]->c;
                double r = diff.norm();
                double m1 = bodies_[i]->density * bodies_[i]->getTemplate().getVolume();
                double m2 = bodies_[j]->density * bodies_[j]->getTemplate().getVolume();
                Fc.segment<3>(3*i) += params_.gravityG*m1*m2/r/r/r * diff;
                Fc.segment<3>(3*j) -= params_.gravityG*m1*m2/r/r/r * diff;
            }
        }
    }
}

void SevenHook::computePenaltyForces(VectorXd &Fc, VectorXd &Ftheta, std::set<Collision> &collisions)
{
    for (auto &collision : collisions) {
        int i = collision.body2, j = collision.body1;
        RigidBodyInstance *Bj = bodies_[j];
        RigidBodyInstance *Bi = bodies_[i];
        Vector3d Vj = Bj->getTemplate().getVerts().row(collision.collidingVertex);
        Matrix3d Ri = VectorMath::rotationMatrix(Bi->theta);
        Matrix3d RiT = Ri.transpose();
        Matrix3d Rj = VectorMath::rotationMatrix(Bj->theta);
        Vector3d Ci = Bi->c;  // q_i
        Vector3d Cj = Bj->c;  // q_j
        Vector3d Vj_q = Rj * Vj + Cj;
        Vector3d Vj_i = RiT * (Vj_q - Ci);
        int teti = collision.collidingTet;
        double D = Bi->getTemplate().distance(Vj_i, teti);
        // std::cout << "D = " << D << "\n";
        Vector3d dD = Bi->getTemplate().Ddistance(teti);
        // V(qi, qj) = k/2 * D(Ri^T * (Rj * Vj + qj - qi))^2;
        // dV(qi, qj, 0i, 0j){qj} = k * D(Ri^T * (Rj * Vj + qj - qi)) * dD(Ri^T * (Rj * Vj + qj - qi)) * Ri^T
        // dRot(0j)Vj{0j} = -Rj * [Vj]x * T(0j)
        // dV(qi, qj, 0i, 0j){0j} = k * D(Ri^T * (Rj * Vj + qj - qi)) * dD(Ri^T * (Rj * Vj + qj - qi)) * Ri^T * dRot(0j)Vj{0j}
        // Vt = Rj * Vj + qj - qi
        // dRotT(0i)Vt{0i} = dRot(-0i)Vt{0i}
        //                 = -Rot(-0i) * [Vt]x * T(-0i) * -1
        //                 = RiT * [Vt]x * T(-0i)
        // dV(qi, qj, 0i, 0j){0i} = k * D(Ri^T * Vt) * dD(Ri^T * Vt) * dRotT(0i)Vt{0i}
        // std::cout << "dD = " << dD << "\n";
        Vector3d Fj = -params_.penaltyStiffness * D * Ri * dD;
        Fc.segment<3>(3*j) += Fj;
        Fc.segment<3>(3*i) -= Fj;

        Matrix3d dRjVj = -1.0 * Rj * VectorMath::crossProductMatrix(Vj) * VectorMath::TMatrix(Bj->theta);
        Vector3d F0j = -params_.penaltyStiffness * D * dRjVj.transpose() * Ri * dD;
        Vector3d Vt = Rj * Vj + Cj - Ci;
        Matrix3d dRiTVt = RiT * VectorMath::crossProductMatrix(Vt) * VectorMath::TMatrix(-Bi->theta);
        Vector3d F0i = -params_.penaltyStiffness * D * dRiTVt.transpose() * dD;
        Ftheta.segment<3>(3*i) += F0i;
        Ftheta.segment<3>(3*j) += F0j;
    }
}

double SevenHook::relv(int i, int j, int k, int tet, double *dist)
{
    RigidBodyInstance *Bj = bodies_[j];
    RigidBodyInstance *Bi = bodies_[i];
    Vector3d Vj = Bj->getTemplate().getVerts().row(k);
    Matrix3d Ri = VectorMath::rotationMatrix(Bi->theta);
    Matrix3d RiT = Ri.transpose();
    Matrix3d Rj = VectorMath::rotationMatrix(Bj->theta);
    Vector3d Ci = Bi->c;  // q_i
    Vector3d Cj = Bj->c;  // q_j
    Vector3d Vj_q = Rj * Vj + Cj;
    Vector3d Vj_i = RiT * (Vj_q - Ci);
    double D = Bi->getTemplate().distance(Vj_i, tet);
    // std::cout << "D = " << D << "\n";
    if (dist != nullptr) {
        *dist = D;
    }

    Vector3d dDj = Ri * Bi->getTemplate().Ddistance(tet);
    // std::cout << "dDj = \n" << dDj << "\n";
    Vector3d dDi = -dDj;
    // std::cout << "dDi = \n" << dDi << "\n";

    // Should I add component for w?
    double v = dDj.dot(Bj->cvel) + dDi.dot(Bi->cvel);
    // std::cout << "Vj = \n" << Bj->cvel << "\n";
    // std::cout << "Vi = \n" << Bi->cvel << "\n";
    
    // std::cout << "v = " << v << "\n";
    return v;
}

Vector3d SevenHook::dg_c_j(int i, int j, int k, int tet)
{
    RigidBodyInstance *Bj = bodies_[j];
    RigidBodyInstance *Bi = bodies_[i];
    Vector3d Vj = Bj->getTemplate().getVerts().row(k);
    Matrix3d Ri = VectorMath::rotationMatrix(Bi->theta);
    Matrix3d RiT = Ri.transpose();
    Matrix3d Rj = VectorMath::rotationMatrix(Bj->theta);
    Vector3d Ci = Bi->c;  // q_i
    Vector3d Cj = Bj->c;  // q_j
    Vector3d Vj_q = Rj * Vj + Cj;
    Vector3d Vj_i = RiT * (Vj_q - Ci);
    // double D = Bi->getTemplate().distance(Vj_i, tet);

    Vector3d dDj = Ri * Bi->getTemplate().Ddistance(tet);
    return dDj;
}

void SevenHook::applyImpulses(std::set<Collision> &collisions)
{
    std::map<int, std::map<int, double> > dist;
    std::map<int, std::map<int, Collision> > pairs;
    for (auto &collision : collisions) {
        int i = collision.body2, j = collision.body1;
        int k = collision.collidingVertex, tet = collision.collidingTet;
        double d = 0;
        double v = relv(i, j, k, tet, &d);
        if (v >= 0) continue;

        int a = (i < j) ? i : j;
        int b = (i < j) ? j : i;

        bool insert = true;
        if (pairs.count(a) > 0 &&
            pairs[a].count(b) > 0) {
            if (dist[a][b] < d) {
                insert = false;
            }
        }

        if (insert) {
            // std::cout << "found collision " << j << "," << i << "\n";
            pairs[a][b] = collision;
            dist[a][b] = d;
        }
    }

    for (auto &p : pairs) {
        auto &cols = p.second;
        for (auto &c : cols) {
            const Collision &collision = c.second;
            // Handle collisionDetection
            int i = collision.body2, j = collision.body1;
            int k = collision.collidingVertex, tet = collision.collidingTet;
            RigidBodyInstance *Bj = bodies_[j];
            RigidBodyInstance *Bi = bodies_[i];
            Vector3d Vj = Bj->getTemplate().getVerts().row(k);
            Matrix3d Ri = VectorMath::rotationMatrix(Bi->theta);
            Matrix3d RiT = Ri.transpose();
            Matrix3d Rj = VectorMath::rotationMatrix(Bj->theta);
            Vector3d Ci = Bi->c;  // q_i
            Vector3d Cj = Bj->c;  // q_j
            Vector3d Vj_q = Rj * Vj + Cj;
            Vector3d Vj_i = RiT * (Vj_q - Ci);
            double Mj = Bj->density * Bj->getTemplate().getVolume();
            double Mi = Bi->density * Bi->getTemplate().getVolume();
            Vector3d dg_j = dg_c_j(i, j, k, tet);
            Vector3d dg_i = -dg_j;
            double v = relv(i, j, k, tet, nullptr);
            double alpha = -1.0 * (1 + params_.CoR) * v / (dg_j.dot(dg_j) / Mj + dg_i.dot(dg_i) / Mi);
            Bj->cvel += alpha * dg_j / Mj;
            Bi->cvel += alpha * dg_i / Mi;
            auto MIj_inv = Bj->getTemplate().getInertiaTensor().inverse() / Bj->density;
            auto MIi_inv = Bi->getTemplate().getInertiaTensor().inverse() / Bi->density;
            Bj->w += MIj_inv * VectorMath::crossProductMatrix(Vj) * Rj.transpose() * alpha * dg_j;
            Bi->w += MIi_inv * VectorMath::crossProductMatrix(Vj_i) * Ri.transpose() * alpha * dg_i;
        }
    }
}

bool SevenHook::simulateOneStep()
{   
    time_ += params_.timeStep;
    int nbodies = (int)bodies_.size();

    std::vector<Vector3d> oldthetas;
    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];
        body.c += params_.timeStep*body.cvel;
        Matrix3d Rhw = VectorMath::rotationMatrix(params_.timeStep*body.w);
        Matrix3d Rtheta = VectorMath::rotationMatrix(body.theta);

        Vector3d oldtheta = body.theta;
        body.theta = VectorMath::axisAngle(Rtheta*Rhw);  
        if (body.theta.dot(oldtheta) < 0 && oldtheta.norm() > M_PI/2.0)
        {
            double oldnorm = oldtheta.norm();
            oldtheta = (oldnorm - 2.0*M_PI)*oldtheta/oldnorm;
        }

        oldthetas.push_back(oldtheta);
    }

    std::set<Collision> collisions;
    collisionDetection(bodies_, collisions);

    Eigen::VectorXd cForce(3 * nbodies);
    Eigen::VectorXd thetaForce(3 * nbodies);
    computeForces(cForce, thetaForce);
    
    // TODO compute and add penalty forces
    if (params_.penaltyEnabled) {
        computePenaltyForces(cForce, thetaForce, collisions);
    }


    // TODO apply collision impulses
    if (params_.impulsesEnabled) {
        applyImpulses(collisions);
    }
    
    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {        
        RigidBodyInstance &body = *bodies_[bodyidx];
        Matrix3d Mi = body.getTemplate().getInertiaTensor();

        body.cvel += params_.timeStep*cForce.segment<3>(3*bodyidx)/body.density/body.getTemplate().getVolume();

        Vector3d newwguess(body.w);
        
        int iter = 0;
        for(iter=0; iter<params_.NewtonMaxIters; iter++)
        {
            Vector3d term1 = (-VectorMath::TMatrix(-params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * newwguess;
            Vector3d term2 = (VectorMath::TMatrix(params_.timeStep*body.w).inverse()*VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * body.w;
            Vector3d term3 = params_.timeStep * thetaForce.segment<3>(3*bodyidx);
            Vector3d fval = term1 + term2 + term3;
            if(fval.norm() / body.density / Mi.trace() <= params_.NewtonTolerance)
                break;

            Matrix3d Df = (-VectorMath::TMatrix(-params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(body.theta)).transpose() * Mi * body.density;

            Vector3d deltaw = Df.inverse() * (-fval);
            newwguess += deltaw;
        }
        // std::cout << "Converged in " << iter << " Newton iterations" << std::endl;
        body.w = newwguess;
    }

    return false;
}

void SevenHook::loadScene()
{
    for (RigidBodyInstance *rbi : bodies_)
        delete rbi;
    for (RigidBodyTemplate *rbt : templates_)
        delete rbt;
    bodies_.clear();
    templates_.clear();

    std::string prefix;
    std::string scenefname = std::string("scenes/") + sceneFile_;
    std::ifstream ifs(scenefname);
    if (!ifs)
    {
        // run from the build directory?
        prefix = "../";
        scenefname = prefix + scenefname;        
        ifs.open(scenefname);
        if(!ifs)
            return;
    }
        

    int nbodies;
    ifs >> nbodies;
    for (int body = 0; body < nbodies; body++)
    {
        std::string meshname;
        ifs >> meshname;
        meshname = prefix + std::string("meshes/") + meshname;
        double scale;
        ifs >> scale;
        RigidBodyTemplate *rbt = new RigidBodyTemplate(meshname, scale);
        double rho;
        ifs >> rho;
        Eigen::Vector3d c, theta, cvel, w;
        for (int i = 0; i < 3; i++)
            ifs >> c[i];
        for (int i = 0; i < 3; i++)
            ifs >> theta[i];
        for (int i = 0; i < 3; i++)
            ifs >> cvel[i];
        for (int i = 0; i < 3; i++)
            ifs >> w[i];
        RigidBodyInstance *rbi = new RigidBodyInstance(*rbt, c, theta, cvel, w, rho);
        templates_.push_back(rbt);
        bodies_.push_back(rbi);
    }
}

