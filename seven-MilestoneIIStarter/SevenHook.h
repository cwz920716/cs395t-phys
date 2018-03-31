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
    // Do not forget to change back
    SevenHook() : PhysicsHook(), sceneFile_("1.scn") {}

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
