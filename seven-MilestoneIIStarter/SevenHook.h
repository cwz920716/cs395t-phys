#include "PhysicsHook.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/ViewerData.h>
#include <deque>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include "SimParameters.h"
#include <set>
#include "CollisionDetection.h"

class RigidBodyTemplate;
class RigidBodyInstance;

class SevenHook : public PhysicsHook
{
public:
    SevenHook() : PhysicsHook(), sceneFile_("box.scn") {}

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);

    virtual void initSimulation();

    virtual void mouseClicked(double x, double y, int button)
    {
    }

    virtual void updateRenderGeometry();

    virtual void tick();

    virtual bool simulateOneStep();

    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
    {
        viewer.data().clear();
        viewer.data().set_mesh(renderQ, renderF);
    }

private:
    void loadScene();
    void computeForces(Eigen::VectorXd &Fc, Eigen::VectorXd &Ftheta);
    void computePenaltyForces(Eigen::VectorXd &Fc, Eigen::VectorXd &Ftheta, std::set<Collision> &collisions);    
    void applyImpulses(std::set<Collision> &collisions);
    double relv(int i, int j, int k, int tet, double *dist);
    Eigen::Vector3d dg_c_j(int i, int j, int k, int tet);
    void explode(int body);

    double time_;
    SimParameters params_;
    std::string sceneFile_;

    std::vector<RigidBodyTemplate *> templates_;
    std::vector<RigidBodyInstance *> bodies_;

    std::mutex destroyMutex_;
    std::vector<int> destroyCommands_;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;

};
