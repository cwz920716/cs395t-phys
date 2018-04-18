#include "FluidHook.h"
#include <limits>
#include <map>
#include <igl/unproject_onto_mesh.h>
#include <algorithm>

FluidHook::FluidHook() : PhysicsHook() 
{
    clickedVertex = -1;
    released = false;

    dt = 1e-3;
    constraintIters = 5;

    gravityEnabled = true;
    gravityG = -9.8;
}

void FluidHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
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
    }    
}

bool FluidHook::mouseClicked(igl::opengl::glfw::Viewer &viewer, int button)
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
        // std::cout << "click at [" << me.pos << "]\n";    
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

bool FluidHook::mouseReleased(igl::opengl::glfw::Viewer &viewer, int button)
{
    MouseEvent me;
    me.type = MouseEvent::ME_RELEASED;
    mouseMutex.lock();
    mouseEvents.push_back(me);
    mouseMutex.unlock();
    return true;
}

bool FluidHook::mouseMoved(igl::opengl::glfw::Viewer &viewer, int button)
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
    return true;
}

void FluidHook::tick()
{
    mouseMutex.lock();
    for (MouseEvent me : mouseEvents)
    {
        if (me.type == MouseEvent::ME_CLICKED)
        {            
            curPos = me.pos;
            clickedPos = me.pos;
            clickedVertex = me.vertex;         
            released = false;
        }
        if (me.type == MouseEvent::ME_RELEASED)
        {
            // std::cout << "release at [" << curPos << "]\n";  
            // clickedVertex = -1;
            released = true;
        }
        if (me.type == MouseEvent::ME_DRAGGED)
        {
            curPos = me.pos;
        }
    }
    mouseEvents.clear();
    mouseMutex.unlock();
}


void FluidHook::initSimulation()
{
    N = 100;
    width = 3.0 / N;

    verts_x = N + 3;
    verts_y = N + 3;
    verts = verts_x * verts_y;

    trans_x = -1.0 * N / 2 * width;
    trans_y = -1.0 * N / 2 * width;

    renderQ.resize(verts, 3);

    for (int y = 0; y < verts_y; y++) {
        for (int x = 0; x < verts_x; x++) {
            int vid = y * verts_x + x;
            Vector3d v(x * width + trans_x, y * width + trans_y, 0);
            renderQ.block<1, 3>(vid, 0) = v;
        }
    }

    faces_x = verts_x - 1;
    faces_y = verts_y - 1;
    faces = 2 * faces_x * faces_y;

    renderF.resize(faces, 3);

    for (int y = 0; y < faces_y; y++) {
        for (int x = 0; x < faces_x; x++) {
            int bl = y * verts_x + x;
            int br = bl + 1;
            int tl = bl + verts_x;
            int tr = tl + 1;
            int fid = y * faces_x + x;
            Vector3i f1(bl, br, tr);
            Vector3i f2(bl, tr, tl);
            renderF.block<1, 3>(fid * 2, 0) = f1;
            renderF.block<1, 3>(fid * 2 + 1, 0) = f2;
        }
    }

    dens.resize(faces_y, faces_x);
    dens.setZero();

    renderC.resize(faces, 3);
    renderC.setZero();
    
}

bool FluidHook::simulateOneStep()
{
    // TODO: time integration

    // add sources
    if (clickedVertex > 0 && released) {
        Vector2i start = pos2grid(clickedPos);
        Vector2i stop = pos2grid(curPos);
        int x_min = std::min(start(0), stop(0));
        int x_max = std::max(start(0), stop(0));
        int y_min = std::min(start(1), stop(1));
        int y_max = std::max(start(1), stop(1));
        // std::cout << "add source for [" << x_min << "-" << x_max << "]["
        //                                 << y_min << "-" << y_max << "]\n";
        clickedVertex = -1;
        released = false;

        for (int y = y_min; y <= y_max; y++) {
            for (int x = x_min; x <= x_max; x++) {
                dens(y, x) = 1.0;
            }
        }
    }

    return false;
}
