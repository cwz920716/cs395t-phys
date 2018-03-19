#include "SevenHook.h"
#include "RigidBodyTemplate.h"
#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include "igl/opengl/glfw/imgui/ImGuiHelpers.h"

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
}

bool SevenHook::simulateOneStep()
{   
    time_ += params_.timeStep;

    const double PI = 3.1415926535897;

    // TODO: rigid body dynamics
    for (RigidBodyInstance *rbi : bodies_)
    {
        rbi->c += rbi->cvel * params_.timeStep;
        auto hw = params_.timeStep * rbi->w;
        auto old_theta = rbi->theta;
        auto axis1 = rbi->theta; axis1.normalize();
        rbi->theta = VectorMath::axisAngle(VectorMath::rotationMatrix(rbi->theta)
                                           * VectorMath::rotationMatrix(hw));
        auto axis2 = rbi->theta;  axis2.normalize();
        // how to handle when theta > pi
        if ((axis1 - axis2).norm() > 1.9) {
            // std::cout << "prepare to flip from " << old_theta.norm() << " to " << rbi->theta.norm() << "\n";
            // std::cout << "prepare to flip a* from\n" << axis1 << "\nto\n" << axis2 << "\n";
            // rbi->theta = (2 * PI - rbi->theta.norm()) * axis1;
            // axis2 = rbi->theta;  axis2.normalize();
            // std::cout << "flip from " << old_theta.norm() << " to " << rbi->theta.norm() << "\n";
            // std::cout << "flip a* from\n" << axis1 << "\nto\n" << axis2 << "\n";
        }

        Eigen::Vector3d w_next = rbi->w;
        Eigen::Matrix3d MI = rbi->getTemplate().getInertiaTensor();
        Eigen::Matrix3d Thw_Tinv = VectorMath::TMatrix(hw).inverse().transpose();
        Eigen::Vector3d rhs1 = Thw_Tinv * MI * rbi->w;
        Eigen::Vector3d hw_next = -1 * params_.timeStep * w_next;
        Eigen::Matrix3d df = VectorMath::TMatrix(hw_next).inverse().transpose() * MI;
        Eigen::Vector3d lhs = df * w_next;
        Eigen::Vector3d f = lhs - rhs1;
        int nIters = 0;
        while (f.norm() > params_.NewtonTolerance) {
            Vector3d delta = df.colPivHouseholderQr().solve(-f);
            w_next += delta;
            hw_next = -1 * params_.timeStep * w_next;
            df = VectorMath::TMatrix(hw_next).inverse().transpose() * MI;
            lhs = df * w_next;
            f = lhs - rhs1;

            nIters++;
            // std::cout << "f.norm = " << f.norm() << "\n";
            // std::cout << "+++ " << nIters << " iter +++\n\n\n";
            if (nIters >= params_.NewtonMaxIters) {
                break;
            }
        }
        // std::cout << "w(i)=\n" << rbi->w << "\n";
        // std::cout << "w(i+1)=\n" << w_next << "\n";
        rbi->w = w_next;
        // std::cout << "theta(i+1)=\n" << rbi->theta.norm() << "\n";
    }

    if (!params_.gravityEnabled) {
        return false;
    }

    int nbody = bodies_.size();
    for (int i = 0; i < nbody; i++) {
        for (int j = i + 1; j < nbody; j++) {
            RigidBodyInstance *b1 = bodies_[i];
            RigidBodyInstance *b2 = bodies_[j];
            Eigen::Vector3d diff = b1->c - b2->c;
            double dist = diff.norm();
            double k = params_.gravityG * b1->getMass() * b2->getMass() / (dist * dist * dist);
            Eigen::Vector3d F1 = -1 * k * diff;
            Eigen::Vector3d F2 = -1 * F1;
            b1->cvel += params_.timeStep * F1 / b1->getMass();
            b2->cvel += params_.timeStep * F2 / b2->getMass();
        }
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
