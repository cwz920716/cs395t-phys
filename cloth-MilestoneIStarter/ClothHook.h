#include "PhysicsHook.h"
#include <igl/readOBJ.h>
#include <iostream>

struct MouseEvent
{
    enum METype {
        ME_CLICKED,
        ME_RELEASED,
        ME_DRAGGED
    };

    METype type;
    int vertex;
    Eigen::Vector3d pos;
};

class ClothHook : public PhysicsHook
{
public:
    ClothHook();

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);
    
    virtual void tick();

    virtual void initSimulation();

    virtual void updateRenderGeometry()
    {
        renderQ = Q;
        renderF = F;        
    }

    virtual bool simulateOneStep();
    
    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
    {
        viewer.data().set_mesh(renderQ, renderF);
    }

    virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, int button);
    virtual bool mouseReleased(igl::opengl::glfw::Viewer &viewer,  int button);
    virtual bool mouseMoved(igl::opengl::glfw::Viewer &viewer,  int button);
    
private:
    Eigen::MatrixXd origQ;
    Eigen::MatrixXd Q;
    Eigen::MatrixXd Qdot;
    Eigen::MatrixXi F;

    float dt;
    int constraintIters;

    bool gravityEnabled;
    float gravityG;
    bool pinEnabled;
    float pinWeight;

    bool stretchEnabled;
    float stretchWeight;

    bool bendingEnabled;
    float bendingWeight;

    bool pullingEnabled;
    float pullingWeight;

    std::mutex mouseMutex;
    std::vector<MouseEvent> mouseEvents;
    int clickedVertex; // the currently selected vertex (-1 if no vertex)
    double clickedz;
    Eigen::Vector3d curPos; // the current position of the mouse cursor in 3D

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
};
