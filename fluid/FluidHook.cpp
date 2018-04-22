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
    dens_intensity = 100;
    constraintIters = 20;

    diff = 1e-3;

    gravityEnabled = true;
    gravityG = -9.8;

    dens_prev = &dens_data0;
    dens = &dens_data1;

    vx_prev = &vx_data0;
    vx = &vx_data1;

    vy_prev = &vy_data0;
    vy = &vy_data1;
}

#define SWAP(a, b) \
do {\
    auto tmp = a; \
    a = b; \
    b = tmp; \
} while(0)

void FluidHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputFloat("Timestep", &dt, 0, 0, 3);        
        ImGui::InputFloat("Diffusion", &diff, 0, 0, 3);        
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
        ret = true;
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
    std::cout << "initSimulation called\n";

    L = 3.0;
    N = 100;
    width = L / N;

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

    renderC.resize(faces, 3);
    renderC.setZero();

    dens_data0.resize(faces_y, faces_x);
    dens_data0.setZero();

    dens_data1.resize(faces_y, faces_x);
    dens_data1.setZero();

    vx_data0.resize(faces_y, faces_x);
    vx_data0.setZero();

    vx_data1.resize(faces_y, faces_x);
    vx_data1.setZero();

    vy_data0.resize(faces_y, faces_x);
    vy_data0.setZero();

    vy_data1.resize(faces_y, faces_x);
    vy_data1.setZero();

    for (int y = 1; y <= N; y++) {
        for (int x = 1; x <= N; x++) {
            (*vy)(y, x) = gravityG;
        }
    }

    // char ch;
    // std::cout << "Press ENTER to continue...\n";
    // std::cin.ignore();
    
    // weird init cond.
    // (*dens_prev)(N/2, N/2) = 100;
}

void FluidHook::dens_step() {
    add_dens_source(*dens, *dens_prev);
    // std::cout << "dens = [\n" << *dens << "]\n";
    // std::cout << "Before diff==============\n";
    SWAP(dens, dens_prev);
    diffuse(0, *dens, *dens_prev);
    SWAP(dens, dens_prev);
    advect(0, *dens, *dens_prev, *vx, *vy);
    // std::cout << "After diff==============\n";
    // std::cout << "dens_s = [\n" << dens_s << "]\n";

    // std::cout << "dens = [\n" << dens << "]\n";
}

void FluidHook::get_sources_from_UI(MatrixXd &d, MatrixXd &u, MatrixXd &v) {
    // d = *dens;
    d.setZero();

    if (clickedVertex > 0 && released) {
        Vector2i start = pos2grid(clickedPos);
        Vector2i stop = pos2grid(curPos);
        int x_min = std::min(start(0), stop(0));
        int x_max = std::max(start(0), stop(0));
        int y_min = std::min(start(1), stop(1));
        int y_max = std::max(start(1), stop(1));
        std::cout << "add source for [" << x_min << "-" << x_max << "]["
                                        << y_min << "-" << y_max << "]\n";
        clickedVertex = -1;
        released = false;

        for (int y = y_min; y <= y_max; y++) {
            for (int x = x_min; x <= x_max; x++) {
                d(y, x) = 1;
            }
        }
    }

}

void FluidHook::add_dens_source(MatrixXd &d, MatrixXd &d0) {
    d += d0;
    for (int y = 0; y < faces_y; y++) {
        for (int x = 0; x < faces_x; x++) {
            if (d(y, x) > 1.0) d(y, x) = 1.0;
            if (d(y, x) < 0.0) d(y, x) = 0.0;
        }
    }
}

void FluidHook::set_bnd(int b, MatrixXd &d) {
    for (int i = 1; i <= N; i++) {
        d(i, 0) = (b == 1) ? -d(i, 1) : d(i, 1);
        d(i, N+1) = (b == 1) ? -d(i, N) : d(i, N);
        d(0, i) = (b == 2) ? -d(1, i) : d(1, i);
        d(N+1, i) = (b == 2) ? -d(N, i) : d(N, i);
    }

    d(0, 0) = 0.5 * (d(0, 1) + d(1, 0));
    d(0, N+1) = 0.5 * (d(0, N) + d(1, N+1));
    d(N+1, 0) = 0.5 * (d(N, 0) + d(N+1, 1));
    d(N+1, N+1) = 0.5 * (d(N, N+1) + d(N+1, N));
}

void FluidHook::diffuse(int b, MatrixXd &d, MatrixXd &d0) {
    double a = diff * dt * N * N;
    d = d0;
    // std::cout << "a=" << a << "\n";

    for (int k = 0; k < constraintIters; k++) {
        // std::cout << "d = [\n" << d << "]\n";
        for (int y = 1; y <= N; y++) {
            for (int x = 1; x <= N; x++) {
                double diffused = d(y, x-1) + d(y, x+1) + d(y-1, x) + d(y+1, x);
                d(y, x) = (d0(y, x) + a * diffused) / (1 + 4 * a);
            }
        }

        set_bnd(b, d);
    }
}

void FluidHook::advect(int b, MatrixXd &d, MatrixXd &d0, MatrixXd &u, MatrixXd &v) {
    d.setZero();
    double dt0 = dt / width;
    for (int y = 1; y <= N; y++) {
        for (int x = 1; x <= N; x++) {
            double xo = x - dt0 * u(y, x);
            double yo = y - dt0 * v(y, x);
            xo = std::max(xo, 0.5);
            xo = std::min(xo, N + 0.5);
            yo = std::max(yo, 0.5);
            yo = std::min(yo, N + 0.5);
            int xi = static_cast<int>(xo);
            int xii = xi + 1;
            int yi = static_cast<int>(yo);
            int yii = yi + 1;
            double s1 = xo - xi;
            double s0 = 1 - s1;
            double t1 = yo - yi;
            double t0 = 1 - t1;
            double dx0 = t0 * d0(yi, xi) + t1 * d0(yii, xi);
            double dx1 = t0 * d0(yii, xii) + t1 * d0(yii, xii);
            d(y, x) = s0 * dx0 + s1 * dx1;
        }
    }

    set_bnd(b, d);
}

void FluidHook::add_vel_source(MatrixXd &u, MatrixXd &v, MatrixXd &d0) {
    if (gravityEnabled) {
        for (int y = 1; y <= N; y++) {
            for (int x = 1; x <= N; x++) {
                v(y, x) += d0(y, x) * dt * gravityG;
            }
        }
    }
}

void FluidHook::project(MatrixXd &u, MatrixXd &v, MatrixXd &p, MatrixXd &div) {
    double h = width;
    for (int y = 1; y <= N; y++) {
        for (int x = 1; x <= N; x++) {
            div(y, x) = -0.5 * h * (u(y, x+1) - u(y, x-1) + v(y+1, x) - v(y-1, x));
            p(y, x) = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);

    p = div;
    for (int k = 0; k < constraintIters; k++) {
        for (int y = 1; y <= N; y++) {
            for (int x = 1; x <= N; x++) {
                p(y, x) = (div(y, x) + p(y-1, x) + p(y+1, x) + p(y, x-1) + p(y, x+1)) / 4;
                
            }
        }

        set_bnd(0, p);
    }

    for (int y = 1; y <= N; y++) {
        for (int x = 1; x <= N; x++) {
            u(y, x) -= 0.5 * (p(y, x+1) - p(y, x-1)) / h;
            v(y, x) -= 0.5 * (p(y+1, x) - p(y+1, x)) / h;
        }
    }

    set_bnd(1, u);
    set_bnd(2, v);
}

void FluidHook::vel_step() {
    // add_vel_source(*vx, *vy, *dens);
    SWAP(vx, vx_prev);
    SWAP(vy, vy_prev);
    diffuse(1, *vx, *vx_prev);
    diffuse(2, *vy, *vy_prev);
    project(*vx, *vy, *vx_prev, *vy_prev);
    SWAP(vx, vx_prev);
    SWAP(vy, vy_prev);
    advect(1, *vx, *vx_prev, *vx_prev, *vy_prev);
    advect(2, *vy, *vy_prev, *vx_prev, *vy_prev);
    project(*vx, *vy, *vx_prev, *vy_prev);
}

bool FluidHook::simulateOneStep()
{
    // TODO: time integration

    get_sources_from_UI(*dens_prev, *vx, *vy);
    // std::cout << "dens_prev=[\n" << *dens_prev << "]\n";
    std::cout << "vy=[\n" << *vy << "]\n";
    vel_step();
    dens_step();
    dens_prev->setZero();

    // char ch;
    // std::cout << "Press ENTER to continue...\n";
    // std::cin.ignore();

    return false;
}
