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

    // explicit Euler for testing
    int n = particles_.size();
    VectorXd q_next(n * 2);
    VectorXd v_next(n * 2);
    VectorXd q = configVector();
    VectorXd v = configVelVector();
    VectorXd q_prev = prevConfigVector();

    // std::cout << "\n==iter begins==\n";

    // std::cout << "G = [\n" << gravity() << "]\n";
    // std::cout << "SpringF = [\n" << springForce(q) << "]\n";

    // std::cout << "F = [\n" << F << "]\n";
    // TODO(wcui): more integrators!
    if (params_.integrator == SimParameters::TimeIntegrator::TI_EXPLICIT_EULER) {
        VectorXd F = gravity() + springForce(q) +
                     viscousDamping(q, q_prev, params_.timeStep) +
                     floorForce(q, q_prev, params_.timeStep);
        // std::cout << "visc = [\n" << viscousDamping(vv) << "]\n";
        // std::cout << "visc = [\n" << viscousDamping(q, q_prev, params_.timeStep) << "]\n";
        // std::cout << "H = [\n" << floorForceHeissan(q, q_prev, params_.timeStep) << "]\n";
        q_next = q + params_.timeStep * v;
        v_next = v + params_.timeStep * massInvMatrix() * F;
    } else if (params_.integrator == SimParameters::TimeIntegrator::TI_VELOCITY_VERLET) {
        q_next = q + params_.timeStep * v;
        VectorXd F = gravity() + springForce(q_next) +
                     viscousDamping(q_next, q, params_.timeStep) +
                     floorForce(q_next, q, params_.timeStep);
        // std::cout << "visc = [\n" << viscousDamping(vv) << "]\n";
        v_next = v + params_.timeStep * massInvMatrix() * F;
    }

    // populate from q_next/v-next to per node pos/vel
    for (int i = 0; i < n; i++) {
        if (particles_[i].fixed)
            continue;

        int i1 = 2 * i, i2 = 2 * i + 1;
        particles_[i].prevpos = particles_[i].pos;
        particles_[i].pos[0] = q_next[i1];
        particles_[i].pos[1] = q_next[i2];
        particles_[i].vel[0] = v_next[i1];
        particles_[i].vel[1] = v_next[i2];
    }

    // TODO(wcui): Handle snapping and sawing
    removeObjects();

    // std::cout << "q = [\n" << q << "];\n\n";
    // std::cout << "q' = [\n" << q_next << "];\n\n";

    // std::cout << "\n==iter ends==\n";
    time_ += params_.timeStep;
    return false;
}

static double shortestP2L(Vector2d p, Vector2d q1, Vector2d q2)
{
    Vector2d ray = q2 - q1;
    if (ray.norm() < 1.0e9) {
        return (p - q1).norm();
    }

    double t0 = (p - q1).dot(ray) / ray.dot(ray);
    Vector2d d = p - (q1 + t0 * ray);
    return d.norm();
}

static double shortestP2S(Vector2d p, Vector2d q1, Vector2d q2)
{
    Vector2d ray = q2 - q1;
    if (ray.norm() < 1.0e-9) {
        return (p - q1).norm();
    }

    double t0 = (p - q1).dot(ray) / ray.dot(ray);
    std::cout << "p=" << p << ", q1=" << q1 << ", q2=" << q2 << "\n";
    std::cout << "t0=" << t0 << "\n";
    if (t0 <= 0) {
        return (p - q1).norm();
    } else if (t0 >= 1) {
        return (p - q2).norm();
    } else {
        Vector2d d = p - (q1 + t0 * ray);
        return d.norm();
    }
}

void GooHook::removeObjects()
{
    std::vector<Particle, Eigen::aligned_allocator<Particle> > live_p;
    std::vector<Connector *> live_s;
    std::vector<Connector *> dead_s;
    std::map<int, int> id_map;

    // VectorXd q = configVector();
    for (int i = 0; i < particles_.size(); i++) {
        // TODO(wcui): remove particles which are too far away
        if (particles_[i].pos.norm() >= 1.42) {
            std::cout << "p left screen\n";
            id_map[i] = -1;
            continue;
        }

        Vector2d pos = particles_[i].pos;
        for (auto saw : saws_) {
            Vector2d diff = pos -saw.pos;
            if (diff.norm() <= saw.radius) {
                std::cout << "particle sawed\n";
                id_map[i] = -1;
                break;
            } 
        }

        if (id_map[i] != -1) {
            int new_id = live_p.size();
            id_map[i] = new_id;
            live_p.push_back(particles_[i]);
        }
    }

    for (auto c : connectors_) {
        auto s = reinterpret_cast<Spring *>(c);
        if (s == nullptr) continue;

        // if one of spring ends is dead, spring is dead
        if (id_map[s->p1] == -1 || id_map[s->p2] == -1) {
            std::cout << "endpoint dead\n";
            dead_s.push_back(s);
            continue;
        }

        // std::cout << "see me?\n";

        // if spring is too long, remove it
        Vector2d q1 = particles_[s->p1].pos;
        Vector2d q2 = particles_[s->p2].pos;
        Vector2d s_vec = q2 - q1;
        double e = s_vec.norm() / s->restlen - 1;
        // std::cout << "e=" << e << "\n";
        if (s->canSnap && e >= params_.maxSpringStrain) {
            std::cout << "snapped\n";
            dead_s.push_back(s);
            continue;
        }

        // if spring intercepts saw, remove it
        bool hit_by_saw = false;
        for (auto saw : saws_) {
            double dist = shortestP2S(saw.pos, q1, q2);
            if (dist <= saw.radius) {
                std::cout << "spring sawed\n";
                hit_by_saw = true;
                dead_s.push_back(s);
                break;
            }
        }

        if (!hit_by_saw) {
            s->p1 = id_map[s->p1];
            s->p2 = id_map[s->p2];
            live_s.push_back(s);
        }
    }

    for (auto c : dead_s) {
        delete c;
    }
    particles_ = live_p;
    connectors_ = live_s;
}

void GooHook::addParticle(double x, double y)
{
    Eigen::Vector2d newpos(x,y);
    double mass = params_.particleMass;
    if(params_.particleFixed)
        mass = std::numeric_limits<double>::infinity();

    int newid = particles_.size();
    particles_.push_back(Particle(newpos, mass, params_.particleFixed, false));

    // TODO connect particles with springs
    for (int i = 0; i < newid; i++) {
        auto &p = particles_[i];
        Vector2d diff = newpos - p.pos;
        double dist = diff.norm();
        if (dist <= params_.maxSpringDist) {
            std::cout << "L=" << dist << "\n";
            Spring *new_spring = new Spring(i, newid, mass + p.mass,
                                            params_.springStiffness,
                                            dist, true);
            connectors_.push_back(new_spring);
        }
    }

    std::cout << "q = [\n" << configVector() << "];\n\n";

    // std::cout << "M = [\n" << massMatrix() << "];\n\n";
}

void GooHook::addSaw(double x, double y)
{
    saws_.push_back(Saw(Vector2d(x,y), params_.sawRadius));    
}

Eigen::VectorXd GooHook::configVector()
{
    Eigen::VectorXd q(particles_.size() * 2);
    for (int i = 0; i < particles_.size(); i++) {
        q[i * 2] = particles_[i].pos[0];
        q[i * 2 + 1] = particles_[i].pos[1];
    }

    return q;
}

Eigen::VectorXd GooHook::prevConfigVector()
{
    Eigen::VectorXd q(particles_.size() * 2);
    for (int i = 0; i < particles_.size(); i++) {
        q[i * 2] = particles_[i].prevpos[0];
        q[i * 2 + 1] = particles_[i].prevpos[1];
    }

    return q;
}

Eigen::VectorXd GooHook::configVelVector()
{
    Eigen::VectorXd q(particles_.size() * 2);
    for (int i = 0; i < particles_.size(); i++) {
        q[i * 2] = particles_[i].vel[0];
        q[i * 2 + 1] = particles_[i].vel[1];
    }

    return q;
}

Eigen::VectorXd GooHook::gravity()
{
    Eigen::VectorXd G(particles_.size() * 2);
    if (!params_.gravityEnabled) {
        G.setZero();
        return G;
    }

    for (int i = 0; i < particles_.size(); i++) {
        G[i * 2] = 0;
        G[i * 2 + 1] = particles_[i].mass * params_.gravityG;
        // std::cout << particles_[i].mass * params_.gravityG;
    }

    return G;
}

Eigen::MatrixXd GooHook::gravityHeissan()
{
    Eigen::MatrixXd dG(particles_.size() * 2, particles_.size() * 2);
    dG.setZero();

    return dG;
}

SpMat GooHook::selector(int i)
{
    SpMat S(2, particles_.size() * 2);
    std::vector<T> vec;
    vec.push_back(T(0, 2 * i, 1));
    vec.push_back(T(1, 2 * i + 1, 1));
    S.setFromTriplets(vec.begin(), vec.end());
    return S;
}

Eigen::MatrixXd GooHook::springForceHeissan(Eigen::VectorXd q)
{
    Eigen::MatrixXd dF(particles_.size() * 2, particles_.size() * 2);
    // TODO(wcui): compute dF for all springs
    return dF;
}

Eigen::VectorXd GooHook::springForce(Eigen::VectorXd q)
{
    Eigen::VectorXd F(particles_.size() * 2);
    F.setZero();
    if (!params_.springsEnabled) {
        return F;
    }
    // Eigen::VectorXd q = configVector();

    for (auto c : connectors_) {
        auto s = reinterpret_cast<Spring *>(c);
        if (s == nullptr) continue;

        SpMat S1 = selector(s->p1);
        SpMat S2 = selector(s->p2);
        Eigen::Vector2d q1 = S1 * q;
        Eigen::Vector2d q2 = S2 * q;
        Eigen::Vector2d diff = q1 - q2;
        double k = s->stiffness / s->restlen;
        // std::cout << "k=" << k << ", L=" << s->restlen << ", |q|=" << diff.norm() << "\n";
        // std::cout << "q=" << diff << "\n";
        
        if (fabs(diff.norm()) <= 1.0e-9) continue;

        double kFactor = k * (diff.norm() - s->restlen) / diff.norm();
        Eigen::VectorXd F_s = kFactor  * (S2.transpose() - S1.transpose()) * diff;
        F += F_s;
    }

    return F;
}


Eigen::VectorXd GooHook::viscousDamping(Eigen::VectorXd q,
                                        Eigen::VectorXd q_prev, double h)
{
    VectorXd v = (q - q_prev) / h;
    Eigen::VectorXd F(particles_.size() * 2);
    F.setZero();
    if (!params_.dampingEnabled) {
        return F;
    }

    for (auto c : connectors_) {
        auto s = reinterpret_cast<Spring *>(c);
        if (s == nullptr) continue;

        Eigen::VectorXd f(particles_.size() * 2);
        f.setZero();

        int p1x = s->p1 * 2;
        int p1y = s->p1 * 2 + 1;
        int p2x = s->p2 * 2;
        int p2y = s->p2 * 2 + 1;

        double kDamp = params_.dampingStiffness;
        f(p1x) = kDamp * (v(p2x) - v(p1x));
        f(p1y) = kDamp * (v(p2y) - v(p1y));
        f(p2x) = kDamp * (v(p1x) - v(p2x));
        f(p2y) = kDamp * (v(p1y) - v(p2y));

        F += f;
    }

    return F;
}

static Eigen::MatrixXd fullSelector(int n, int i, int o)
{
    Eigen::MatrixXd S(2 * n, 2 * n);
    S.setZero();
    S(2 * o, 2 * i) = 1;
    S(2 * o + 1, 2 * i + 1) = 1;
    return S;
}

Eigen::MatrixXd GooHook::viscousDampingHeissan(Eigen::VectorXd q,
                                               Eigen::VectorXd q_prev, double h)
{
    int n = particles_.size();
    Eigen::MatrixXd dF(n * 2, n * 2);
    dF.setZero();

    for (auto c : connectors_) {
        auto s = reinterpret_cast<Spring *>(c);
        if (s == nullptr) continue;

        double k = params_.dampingStiffness / h;

        dF += k * (fullSelector(n, s->p2, s->p1) - fullSelector(n, s->p1, s->p1));
        dF += k * (fullSelector(n, s->p1, s->p2) - fullSelector(n, s->p2, s->p2));
    }

    return dF;
}

Eigen::VectorXd GooHook::floorForce(Eigen::VectorXd q,
                                    Eigen::VectorXd q_prev, double h)
{
    VectorXd v = (q - q_prev) / h;
    Eigen::VectorXd F(particles_.size() * 2);
    if (!params_.floorEnabled) {
        F.setZero();
        return F;
    }

    // potential part
    for (int i = 0; i < particles_.size(); i++) {
        F[i * 2] = 0;
        if (q[i * 2 + 1] < -0.5) {
            F[i * 2 + 1] = -floorImpulse * particles_[i].mass * params_.gravityG;
        } else {
            F[i * 2 + 1] = 0;
        }
    }

    // non-conservative part, the drag force
    double kDrag = -5.0 * floorDrag;
    // F += kDrag * massMatrix() * v;
    for (int i = 0; i < particles_.size(); i++) {
        if (q[i * 2 + 1] < -0.5) {
            F[i * 2 + 0] += kDrag  * particles_[i].mass * v[i * 2 + 0];
            F[i * 2 + 1] += kDrag  * particles_[i].mass * v[i * 2 + 1];
        }
    }

    return F;
}

Eigen::MatrixXd GooHook::floorForceHeissan(Eigen::VectorXd q,
                                           Eigen::VectorXd q_prev, double h)
{
    int n = particles_.size();
    Eigen::MatrixXd dF(n * 2, n * 2);
    Eigen::MatrixXd I(n * 2, n * 2);
    dF.setZero();
    I.setIdentity();

    // TODO: fix it!
    // double k = -floorDrag / h;
    // dF += k * I;

    return dF;
}

SpMat GooHook::massMatrix()
{
    int n = particles_.size();
    SpMat M(n * 2, n * 2);
    std::vector<T> mass_vec;
    for (int i = 0; i < n; i++) {
        int i1 = i * 2, i2 = i * 2 + 1;
        mass_vec.push_back(T(i1, i1, particles_[i].mass));
        mass_vec.push_back(T(i2, i2, particles_[i].mass));
    }

    M.setFromTriplets(mass_vec.begin(), mass_vec.end());
    return M;
}

SpMat GooHook::massInvMatrix()
{
    int n = particles_.size();
    SpMat M_inv(n * 2, n * 2);
    std::vector<T> mass_inv_vec;
    for (int i = 0; i < n; i++) {
        int i1 = i * 2, i2 = i * 2 + 1;
        double mass_inv = 1.0 / particles_[i].mass;
        mass_inv_vec.push_back(T(i1, i1, mass_inv));
        mass_inv_vec.push_back(T(i2, i2, mass_inv));
    }

    M_inv.setFromTriplets(mass_inv_vec.begin(), mass_inv_vec.end());
    return M_inv;
}
