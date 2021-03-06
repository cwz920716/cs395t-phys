#include "PhysicsHook.h"
#include <igl/readOBJ.h>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <Eigen/SVD>

using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

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
    Eigen::MatrixXd C0;
    Eigen::MatrixXi Quads;
    Eigen::MatrixXd QC0;

    float dt;
    int constraintIters;

    bool gravityEnabled;
    float gravityG;
    bool pinEnabled;
    float pinWeight;

    int pinTL, pinTR;
    Eigen::Vector3d xTL, xTR;
    SpMat selector(int i);
    SpMat selectorT(int i);
    void applyStretch(Eigen::VectorXd &Q);
    void applyBending(Eigen::VectorXd &Q);

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
