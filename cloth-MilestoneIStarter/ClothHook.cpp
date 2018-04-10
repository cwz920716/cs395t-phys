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
    double xmin = Q(0, 0), xmax = Q(0, 0), ymin = Q(0, 1), ymax = Q(0, 1);
    pinTL = pinTR = 0;
    for (int i = 0; i < Q.rows(); i++) {
        if (Q(i, 0) <= xmin && Q(i, 1) >= ymax) {
            pinTL = i;
        }
        if (Q(i, 0) >= xmax && Q(i, 1) >= ymax) {
            pinTR = i;
        }
        if (Q(i, 0) <= xmin) {
            xmin = Q(i, 0);
        }
        if (Q(i, 0) >= xmax) {
            xmax = Q(i, 0);
        }
        if (Q(i, 1) <= ymin) {
            ymin = Q(i, 1);
        }
        if (Q(i, 1) >= ymax) {
            ymax = Q(i, 1);
        }
    }

    std::cout << "[" << xmin << "~~" << xmax << "]["
              << ymin << "~~" << ymax << "]\n";
    xTL = Q.block<1,3>(pinTL, 0);
    xTR = Q.block<1,3>(pinTR, 0);
    std::cout << "pinTL= [\n" << xTL << "]\n @@" << pinTL << "\n";
    std::cout << "pinTR= [\n" << xTR << "]\n @@" << pinTR << "\n";
    std::cout << "F=[" << F.rows() << ", " << F.cols() << "]\n";

    C0.resize(F.rows(), 3);
    C0.setZero();
    for (int i = 0; i < F.rows(); i++) {
        int a = F(i, 0), b = F(i, 1), c = F(i, 2);
        Eigen::Vector3d A, B, C;
        for (int j = 0; j < 3; j++) {
            A(j) = origQ(a, j);
            B(j) = origQ(b, j);
            C(j) = origQ(c, j);
        }

        // std::cout << "A=" << A << "\n";
        // std::cout << "B=" << B << "\n";
        // std::cout << "C=" << C << "\n";

        Eigen::Vector3d center = (A + B + C) / 3.0;
        for (int j = 0; j < 3; j++) {
            C0(i, j) = center(j);
        }
        // std::cout << "C0[i] = " << C0.block<1, 3>(i, 0) << "\n";
    }

    
}

SpMat ClothHook::selector(int i)
{
    SpMat S(3, Q.rows() * 3);
    std::vector<T> vec;
    vec.push_back(T(0, 3 * i, 1));
    vec.push_back(T(1, 3 * i + 1, 1));
    vec.push_back(T(2, 3 * i + 2, 1));
    S.setFromTriplets(vec.begin(), vec.end());
    return S;
}

SpMat ClothHook::selectorT(int i)
{
    SpMat S(3, Q.rows() * 3);
    std::vector<T> vec;
    vec.push_back(T(0, 3 * i, 1));
    vec.push_back(T(1, 3 * i + 1, 1));
    vec.push_back(T(2, 3 * i + 2, 1));
    S.setFromTriplets(vec.begin(), vec.end());
    return S.transpose();
}

void ClothHook::applyStretch(Eigen::VectorXd &q) {
    std::map<int, bool> projected;
    for (int i = 0; i < F.rows(); i++) {
        int x0 = F(i, 0), x1 = F(i, 1), x2 = F(i, 2);
        Eigen::Vector3d X0, X1, X2, center0, O0, O1, O2;
        for (int j = 0; j < 3; j++) {
            X0(j) = q(3 * x0 + j);
            X1(j) = q(3 * x1 + j);
            X2(j) = q(3 * x2 + j);
            center0(j) = C0(i, j);
            O0(j) = origQ(x0, j);
            O1(j) = origQ(x1, j);
            O2(j) = origQ(x2, j);
        }
        // std::cout << "X0 = [" << X0 << "]\n";
        // std::cout << "X1 = [" << X1 << "]\n";
        // std::cout << "X2 = [" << X2 << "]\n";
        Eigen::Vector3d center = (X0 + X1 + X2) / 3.0;
        // std::cout << "c0 = " << center0 << "\n";
        // std::cout << "c0' = " << center << "\n";
        Eigen::Vector3d T0 = center - center0;
        // std::cout << "Trans = [" << T0.norm() << "]\n";

        Eigen::Matrix3d A, B;
        A.block<3, 1>(0, 0) = X0 - center;
        A.block<3, 1>(0, 1) = X1 - center;
        A.block<3, 1>(0, 2) = X2 - center;
        B.block<3, 1>(0, 0) = O0 - center0;
        B.block<3, 1>(0, 1) = O1 - center0;
        B.block<3, 1>(0, 2) = O2 - center0;
        // std::cout << "A = [" << A << "]\n";
        // std::cout << "B = [" << B << "]\n";

        auto M = A * B.transpose();
        JacobiSVD<MatrixXd> svd(M, ComputeFullU | ComputeFullV);
        Matrix3d R = svd.matrixU() * svd.matrixV().transpose();
        // std::cout << R << "\n for triangle " << i << "\n";
        auto Px0 = R * (O0 - center0) + center;
        auto Px1 = R * (O1 - center0) + center;
        auto Px2 = R * (O2 - center0) + center;

        // projected[x0] = projected[x1] = projected[x2] = true;

        q.segment<3>(3 * x0) = stretchWeight * Px0 + (1 - stretchWeight)  * q.segment<3>(3 * x0);
        q.segment<3>(3 * x1) = stretchWeight * Px1 + (1 - stretchWeight)  * q.segment<3>(3 * x1);
        q.segment<3>(3 * x2) = stretchWeight * Px2 + (1 - stretchWeight)  * q.segment<3>(3 * x2);
    }

    // q = stretchWeight * Pq + (1 - stretchWeight)  * q;
}

bool ClothHook::simulateOneStep()
{
    // TODO: time integration
    // std::cout << "|v|=[" << Qdot.norm() << "]\n";

    auto Qold = Q;
    Q = Qold + dt * Qdot;

    Eigen::VectorXd Qnext(Q.rows()*Q.cols());
    Qnext.setZero();
    for (int i = 0; i < Q.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            Qnext(3 * i + j) = Q(i, j);
        }
    }
    auto coq = Qnext;

    for (int iter = 0; iter < constraintIters; iter++) {
        // std::cout << "---" << iter << "---\n";

        if (pinEnabled) {
            // std::cout << pinTL << "\n";
            // std::cout << "Si = \n" << selector(pinTL) << "\n";
            Eigen::VectorXd Qproj = Qnext;// - selectorT(pinTL) * selector(pinTL) * Qnext + selectorT(pinTL) * xTL;
            Qproj.segment<3>(3 * pinTL) = xTL;
            Qnext = pinWeight * Qproj + (1 - pinWeight) * Qnext;
            // std::cout << Qnext.segment<3>(3 * pinTL);
            // std::cout << xTL;
            Qproj = Qnext;
            Qproj.segment<3>(3 * pinTR) = xTR;
            Qnext = pinWeight * Qproj + (1 - pinWeight) * Qnext;
        }
        // std::cout << "|dq|=[" << (Qnext - coq).norm() << "]\n";

        if (stretchEnabled) {
            applyStretch(Qnext);
        }
        // std::cout << "|dq|=[" << (Qnext - coq).norm() << "]\n";

        if (pullingEnabled && clickedVertex >= 0 && clickedVertex < Q.rows()) {
            Eigen::VectorXd Qproj = Qnext;
            Qproj.segment<3>(3 * clickedVertex) = curPos;
            Qnext = pullingWeight * Qproj + (1 - pullingWeight) * Qnext;
        }
    }

    for (int i = 0; i < Q.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            Q(i, j) = Qnext(3 * i + j);
        }
    }

    Qdot = (Q - Qold) / dt;
    int Down = 1;
    if (gravityEnabled) {
        for (int i = 0; i < Qdot.rows(); i++) {
            Qdot(i, Down) += dt * gravityG;
        }
    }

    return false;
}
