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
    // TODO apply collision impulses
    
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
        std::cout << "Converged in " << iter << " Newton iterations" << std::endl;
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

