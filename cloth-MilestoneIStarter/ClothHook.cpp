#include "ClothHook.h"
#include <limits>
#include <map>
#include <igl/unproject_onto_mesh.h>

ClothHook::ClothHook() : PhysicsHook() 
{
    clickedVertex = -1;

    dt = 1e-3;
    constraintIters = 5;

    gravityEnabled = true;
    gravityG = -9.8;

    pinEnabled = true;
    pinWeight = 1.0;

    stretchEnabled = true;
    stretchWeight = 0.5;

    bendingEnabled = true;
    bendingWeight = 0.5;

    pullingEnabled = true;
    pullingWeight = 0.5;
}

void ClothHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputFloat("Timestep", &dt, 0, 0, 3);        
        ImGui::InputInt("Constraint Iters", &constraintIters);
    }
    if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Gravity Enabled", &gravityEnabled);
        ImGui::InputFloat("Gravity G", &gravityG, 0, 0, 3);        
        ImGui::Checkbox("Pins Enabled", &pinEnabled);
        ImGui::InputFloat("Pin Weight", &pinWeight, 0, 0, 3);        
        ImGui::Checkbox("Stretching Enabled", &stretchEnabled);
        ImGui::InputFloat("Stretching Weight", &stretchWeight, 0, 0, 3);    
        ImGui::Checkbox("Bending Enabled", &bendingEnabled);
        ImGui::InputFloat("Bending Weight", &bendingWeight, 0, 0, 3);   
        ImGui::Checkbox("Pulling Enabled", &pullingEnabled);
        ImGui::InputFloat("Pulling Weight", &pullingWeight, 0, 0, 3);   
    }    
}

bool ClothHook::mouseClicked(igl::opengl::glfw::Viewer &viewer, int button)
{
    if(button != 0)
        return false;
    render_mutex.lock();
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    int fid;
    Eigen::Vector3f bc;
    MouseEvent me;
    bool ret = true;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
        viewer.core.proj, viewer.core.viewport, renderQ, renderF, fid, bc))
    {
        int bestvert = -1;
        double bestcoord = 2.0;
        for (int j = 0; j < 3; j++)
        {
            if (bc[j] < bestcoord)
            {
                bestcoord = bc[j];
                bestvert = j;
            }
        }
        me.type = MouseEvent::ME_CLICKED;
        me.vertex = renderF(fid, bestvert);        
        
        Eigen::Vector3f proj;
        Eigen::Vector3f pt;
        for (int i = 0; i < 3; i++)
            pt[i] = renderQ(me.vertex, i);
        Eigen::Matrix4f modelview = viewer.core.view * viewer.core.model;
        proj = igl::project(pt, modelview,
            viewer.core.proj, viewer.core.viewport);

        clickedz = proj[2];
        Eigen::Vector3f pos;
        pos[0] = float(x);
        pos[1] = float(y);
        pos[2] = float(clickedz);
        Eigen::Vector3f unproj = igl::unproject(pos, modelview,
            viewer.core.proj, viewer.core.viewport);
        for (int i = 0; i < 3; i++)
            me.pos[i] = unproj[i];
        ret = true;
    }
    else
    {
        me.type = MouseEvent::ME_RELEASED;
        ret = false;
    }       
    render_mutex.unlock();

    mouseMutex.lock();
    mouseEvents.push_back(me);
    mouseMutex.unlock();
    return ret;
}

bool ClothHook::mouseReleased(igl::opengl::glfw::Viewer &viewer, int button)
{
    MouseEvent me;
    me.type = MouseEvent::ME_RELEASED;
    mouseMutex.lock();
    mouseEvents.push_back(me);
    mouseMutex.unlock();
    return false;
}

bool ClothHook::mouseMoved(igl::opengl::glfw::Viewer &viewer, int button)
{
    MouseEvent me;
    me.type = MouseEvent::ME_DRAGGED;    
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    Eigen::Vector3d pos(x, y, clickedz);
    igl::unproject(pos, viewer.core.view * viewer.core.model,
        viewer.core.proj, viewer.core.viewport, me.pos);
    mouseMutex.lock();
    mouseEvents.push_back(me);
    mouseMutex.unlock();
    return false;
}

void ClothHook::tick()
{
    mouseMutex.lock();
    for (MouseEvent me : mouseEvents)
    {
        if (me.type == MouseEvent::ME_CLICKED)
        {            
            curPos = me.pos;
            clickedVertex = me.vertex;         
        }
        if (me.type == MouseEvent::ME_RELEASED)
        {
            clickedVertex = -1;
        }
        if (me.type == MouseEvent::ME_DRAGGED)
        {
            curPos = me.pos;
        }
    }
    mouseEvents.clear();
    mouseMutex.unlock();
}


void ClothHook::initSimulation()
{
    if(!igl::readOBJ("meshes/rect-coarse.obj", origQ, F))
        if (!igl::readOBJ("../meshes/rect-coarse.obj", origQ, F))
        {
            std::cerr << "Couldn't read mesh file" << std::endl;
            exit(-1);
        }
    //mesh is tiny for some reason
    origQ *= 50;
    Q = origQ;  
    Qdot.resize(Q.rows(), 3);
    Qdot.setZero();

    // TODO: precompute things like the top corners of the cloth, the bending hinges, etc.
}

bool ClothHook::simulateOneStep()
{
    // TODO: time integration
    return false;
}
