#include "PhysicsHook.h"
#include <igl/readOBJ.h>
#include <iostream>

using namespace Eigen;

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

class FluidHook : public PhysicsHook
{
public:
    FluidHook();

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);
    
    virtual void tick();

    virtual void initSimulation();

    virtual void updateRenderGeometry()
    {
        if (dens == nullptr) {
            std::cout << "Error: dens should never be nullptr!\n";
            exit(0);
        }

        for (int y = 0; y < faces_y; y++) {
            for (int x = 0; x < faces_x; x++) {
                int fid = y * faces_x + x;
                Vector3d c(0, 0, (*dens)(y, x));
                renderC.block<1, 3>(fid * 2, 0) = c;
                renderC.block<1, 3>(fid * 2 + 1, 0) = c;
            }
        }
    }

    virtual bool simulateOneStep();
    
    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
    {
        viewer.data().set_mesh(renderQ, renderF);
        viewer.data().set_colors(renderC);
    }

    virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, int button);
    virtual bool mouseReleased(igl::opengl::glfw::Viewer &viewer,  int button);
    virtual bool mouseMoved(igl::opengl::glfw::Viewer &viewer,  int button);
    
private:

    Vector2i pos2grid(Vector3d pos) {
        int x = (pos(0) - trans_x) / width;
        if (x <= 1) x = 1;
        if (x >= N) x = N;
        int y = (pos(1) - trans_y) / width;
        if (y <= 1) y = 1;
        if (y >= N) y = N;
        return Vector2i(x, y); 
    }

    void get_sources_from_UI(MatrixXd &d);
    void dens_step();
    void diffuse(MatrixXd &d, MatrixXd &d0);
    void add_source(MatrixXd &d, MatrixXd &d0);
    void set_bnd(int b, MatrixXd &d);

    float dt;
    float dens_intensity;
    float diff;
    int constraintIters;

    bool gravityEnabled;
    float gravityG;

    int N;
    double width;

    std::mutex mouseMutex;
    std::vector<MouseEvent> mouseEvents;
    int clickedVertex; // the currently selected vertex (-1 if no vertex)
    bool released;
    double clickedz;
    Eigen::Vector3d curPos, clickedPos; // the current position of the mouse cursor in 3D

    Eigen::MatrixXd dens_data0, dens_data1;
    MatrixXd *dens, *dens_prev;

    int verts, verts_x, verts_y;
    int faces, faces_x, faces_y;
    double trans_x, trans_y;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd renderC;
};
