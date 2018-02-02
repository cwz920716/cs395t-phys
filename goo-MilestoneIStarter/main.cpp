#include <igl/viewer/Viewer.h>
#include <thread>
#include "PhysicsHook.h"
#include "GooHook.h"
#include <igl/unproject.h>

static PhysicsHook *hook = NULL;

void toggleSimulation()
{
    if (!hook)
        return;

    if (hook->isPaused())
        hook->run();
    else
        hook->pause();
}

void resetSimulation()
{
    if (!hook)
        return;

    hook->reset();
}

bool drawCallback(igl::viewer::Viewer &viewer)
{
    if (!hook)
        return false;

    hook->render(viewer);
    return false;
}

bool keyCallback(igl::viewer::Viewer& viewer, unsigned int key, int modifiers)
{
    if (key == ' ')
    {
        toggleSimulation();
        return true;
    }
    return false;
}

bool mouseCallback(igl::viewer::Viewer& viewer, int button, int modifier)
{
    Eigen::Vector3f pos(viewer.down_mouse_x, viewer.down_mouse_y, 0);
    Eigen::Matrix4f model = viewer.core.view*viewer.core.model;
    Eigen::Vector3f unproj = igl::unproject(pos, model, viewer.core.proj, viewer.core.viewport);
    hook->mouseClicked(unproj[0], -unproj[1], button);
    return true;
}

bool mouseScroll(igl::viewer::Viewer& viewer, float delta)
{
    return true;
}

bool initGUI(igl::viewer::Viewer &viewer)
{
    viewer.ngui->window()->setVisible(false);
    viewer.ngui->addWindow(Eigen::Vector2i(0, 0), "Simulation Control");
    viewer.ngui->addButton("Run/Pause Sim", toggleSimulation);
    viewer.ngui->addButton("Reset Sim", resetSimulation);
    hook->initGUI(viewer);
    viewer.screen->performLayout();
    return false;
}

int main(int argc, char *argv[])
{
  igl::viewer::Viewer viewer;

  hook = new GooHook();
  hook->reset();
  viewer.core.orthographic = true;
  viewer.core.camera_zoom = 4.0;
  viewer.core.show_lines = false;
  viewer.data.set_face_based(false);
  viewer.core.is_animating = true;
  viewer.callback_key_pressed = keyCallback;
  viewer.callback_pre_draw = drawCallback;
  viewer.callback_mouse_down = mouseCallback;
  viewer.callback_mouse_scroll = mouseScroll;
  viewer.callback_init = initGUI;
  viewer.launch();
}
