#include "GooHook.h"

using namespace Eigen;

void GooHook::initGUI(igl::viewer::Viewer &viewer)
{
    viewer.ngui->addGroup("UI Options");
    viewer.ngui->addVariable("Click Adds", params_.clickMode)->setItems({ "Particles", "Saws" });
    viewer.ngui->addGroup("Simulation Options");
    viewer.ngui->addVariable("Timestep", params_.timeStep);
    viewer.ngui->addVariable("Integrator", params_.integrator)->setItems({ "Explicit Euler", "Implicit Euler", "Implicit Midpoint", "Velocity Verlet" });
    viewer.ngui->addVariable("Newton Tolerance", params_.NewtonTolerance);
    viewer.ngui->addVariable("Newton Max Iters", params_.NewtonMaxIters);
    viewer.ngui->addGroup("Forces");
    viewer.ngui->addVariable("Gravity Enabled", params_.gravityEnabled);
    viewer.ngui->addVariable("  Gravity g", params_.gravityG);
    viewer.ngui->addVariable("Springs Enabled", params_.springsEnabled);
    viewer.ngui->addVariable("  Max Strain", params_.maxSpringStrain);
    viewer.ngui->addVariable("Damping Enabled", params_.dampingEnabled);
    viewer.ngui->addVariable("  Viscosity", params_.dampingStiffness);
    viewer.ngui->addVariable("Floor Enabled", params_.floorEnabled);
    
    viewer.ngui->addWindow(Eigen::Vector2i(1000, 0), "New Objects");
    viewer.ngui->addGroup("New Particles");
    viewer.ngui->addVariable("Is Fixed", params_.particleFixed);
    viewer.ngui->addVariable("Mass", params_.particleMass);
    viewer.ngui->addGroup("New Saws");
    viewer.ngui->addVariable("Radius", params_.sawRadius);
    viewer.ngui->addGroup("New Springs");
    viewer.ngui->addVariable("Max Spring Dist", params_.maxSpringDist);
    viewer.ngui->addVariable("Base Stiffness", params_.springStiffness);    
}

void GooHook::updateRenderGeometry()
{
    double baseradius = 0.02;
    double pulsefactor = 0.1;
    double pulsespeed = 50.0;

    int sawteeth = 20;
    double sawdepth = 0.1;
    double sawangspeed = 10.0;

    double baselinewidth = 0.005;

    int numcirclewedges = 20;

    // this is terrible. But, easiest to get up and running

    std::vector<Eigen::Vector3d> verts;
    std::vector<Eigen::Vector3d> vertexColors;
    std::vector<Eigen::Vector3i> faces;

    int idx = 0;

    double eps = 1e-4;


    if(params_.floorEnabled)
    {
        for (int i = 0; i < 6; i++)
        {
            vertexColors.push_back(Eigen::Vector3d(0.3, 1.0, 0.3));
        }

        verts.push_back(Eigen::Vector3d(-1, -0.5, eps));
        verts.push_back(Eigen::Vector3d(1, -0.5, eps));
        verts.push_back(Eigen::Vector3d(-1, -1, eps));

        faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));

        verts.push_back(Eigen::Vector3d(-1, -1, eps));
        verts.push_back(Eigen::Vector3d(1, -0.5, eps));
        verts.push_back(Eigen::Vector3d(1, -1, eps));
        faces.push_back(Eigen::Vector3i(idx + 3, idx + 4, idx + 5));
        idx += 6;
    }

    
    for (std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
    {
        Eigen::Vector3d color;
        if ((*it)->associatedBendingStencils.empty())
            color << 0.0, 0.0, 1.0;
        else
            color << 0.75, 0.5, 0.75;
        Vector2d sourcepos = particles_[(*it)->p1].pos;
        Vector2d destpos = particles_[(*it)->p2].pos;

        Vector2d vec = destpos - sourcepos;
        Vector2d perp(-vec[1], vec[0]);
        perp /= perp.norm();

        double dist = (sourcepos - destpos).norm();

        double width = baselinewidth / (1.0 + 20.0 * dist * dist);

        for (int i = 0; i < 4; i++)
            vertexColors.push_back(color);

        verts.push_back(Eigen::Vector3d(sourcepos[0] + width * perp[0], sourcepos[1] + width * perp[1], -eps));
        verts.push_back(Eigen::Vector3d(sourcepos[0] - width * perp[0], sourcepos[1] - width * perp[1], -eps));
        verts.push_back(Eigen::Vector3d(destpos[0] + width * perp[0], destpos[1] + width * perp[1], -eps));
        verts.push_back(Eigen::Vector3d(destpos[0] - width * perp[0], destpos[1] - width * perp[1], -eps));

        faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));
        faces.push_back(Eigen::Vector3i(idx + 2, idx + 1, idx + 3));
        idx += 4;
    }

    int nparticles = particles_.size();

    for(int i=0; i<nparticles; i++)
    {
        double radius = baseradius*sqrt(particles_[i].mass);
        radius *= (1.0 + pulsefactor*sin(pulsespeed*time_));

        Eigen::Vector3d color(0,0,0);

        if(particles_[i].fixed)
        {
            radius = baseradius;
            color << 1.0, 0, 0;
        }

        for (int j = 0; j < numcirclewedges + 2; j++)
        {
            vertexColors.push_back(color);
        }


        verts.push_back(Eigen::Vector3d(particles_[i].pos[0], particles_[i].pos[1], 0));
        
        const double PI = 3.1415926535898;
        for (int j = 0; j <= numcirclewedges; j++)
        {
            verts.push_back(Eigen::Vector3d(particles_[i].pos[0] + radius * cos(2 * PI*j / numcirclewedges),
                particles_[i].pos[1] + radius * sin(2 * PI*j / numcirclewedges), 0));            
        }

        for (int j = 0; j <= numcirclewedges; j++)
        {
            faces.push_back(Eigen::Vector3i(idx, idx + j + 1, idx + 1 + ((j + 1) % (numcirclewedges + 1))));
        }
        
        idx += numcirclewedges + 2;
    }

    for(std::vector<Saw>::iterator it = saws_.begin(); it != saws_.end(); ++it)
    {
        double outerradius = it->radius;
        double innerradius = (1.0-sawdepth)*outerradius;

        Eigen::Vector3d color(0.5,0.5,0.5);

        int spokes = 2*sawteeth;
        for (int j = 0; j < spokes + 2; j++)
        {
            vertexColors.push_back(color);
        }

        verts.push_back(Eigen::Vector3d(it->pos[0], it->pos[1], 0));
        
        const double PI = 3.1415926535898;
        for (int i = 0; i <= spokes; i++)
        {
            double radius = (i % 2 == 0) ? innerradius : outerradius;
            verts.push_back(Eigen::Vector3d(it->pos[0] + radius * cos(2 * PI*i / spokes + sawangspeed*time_),
                it->pos[1] + radius * sin(2 * PI*i / spokes + sawangspeed*time_), 0));
        }

        for (int j = 0; j <= spokes; j++)
        {
            faces.push_back(Eigen::Vector3i(idx, idx + j + 1, idx + 1 + ((j + 1) % (spokes + 1))));
        }

        idx += spokes + 2;        
    }

    renderQ.resize(verts.size(),3);
    renderC.resize(vertexColors.size(), 3);
    for (int i = 0; i < verts.size(); i++)
    {
        renderQ.row(i) = verts[i];
        renderC.row(i) = vertexColors[i];
    }
    renderF.resize(faces.size(), 3);
    for (int i = 0; i < faces.size(); i++)
        renderF.row(i) = faces[i];
}


void GooHook::initSimulation()
{
    time_ = 0;
    particles_.clear();
    for(std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
        delete *it;
    connectors_.clear();
    saws_.clear();
}

void GooHook::tick()
{
    message_mutex.lock();
    {
        while (!mouseClicks_.empty())
        {
            MouseClick mc = mouseClicks_.front();
            mouseClicks_.pop_front();            
            switch (mc.mode)
            {
            case SimParameters::ClickMode::CM_ADDPARTICLE:
            {
                addParticle(mc.x, mc.y);
                break;
            }
            case SimParameters::ClickMode::CM_ADDSAW:
            {
                addSaw(mc.x, mc.y);
                break;
            }
            }
        }
    }
    message_mutex.unlock();    
}

bool GooHook::simulateOneStep()
{
    // TODO: implement time integration
    time_ += params_.timeStep;
    return false;
}

void GooHook::addParticle(double x, double y)
{
    Vector2d newpos(x,y);
    double mass = params_.particleMass;
    if(params_.particleFixed)
        mass = std::numeric_limits<double>::infinity();

    int newid = particles_.size();
    particles_.push_back(Particle(newpos, mass, params_.particleFixed, false));

    // TODO connect particles with springs
}

void GooHook::addSaw(double x, double y)
{
    saws_.push_back(Saw(Vector2d(x,y), params_.sawRadius));    
}
