#include "GooHook.h"

using namespace Eigen;

void GooHook::initGUI(igl::viewer::Viewer &viewer)
{
    viewer.ngui->addGroup("UI Options");
    viewer.ngui->addVariable("Click Adds", params_.clickMode)->setItems({ "Particles", "Saws" });
    viewer.ngui->addVariable("Connector Type", params_.connectorType)->setItems({ "Springs", "Rigid Rods", "Flexible Rods" });
    viewer.ngui->addGroup("Simulation Options");
    viewer.ngui->addVariable("Constraint Handling", params_.constraintHandling)->setItems({ "Penalty Method", "Step and Project", "Lagrange Multipliers" });
    viewer.ngui->addVariable("Timestep", params_.timeStep);
    viewer.ngui->addVariable("Newton Tolerance", params_.NewtonTolerance);
    viewer.ngui->addVariable("Newton Max Iters", params_.NewtonMaxIters);
    viewer.ngui->addVariable("Penalty Stiffness", params_.penaltyStiffness);
    viewer.ngui->addGroup("Forces");
    viewer.ngui->addVariable("Gravity Enabled", params_.gravityEnabled);
    viewer.ngui->addVariable("  Gravity g", params_.gravityG);
    viewer.ngui->addVariable("Springs Enabled", params_.springsEnabled);
    viewer.ngui->addVariable("  Max Strain", params_.maxSpringStrain);
    viewer.ngui->addVariable("Damping Enabled", params_.dampingEnabled);
    viewer.ngui->addVariable("  Viscosity", params_.dampingStiffness);
    viewer.ngui->addVariable("Floor Enabled", params_.floorEnabled);
    viewer.ngui->addVariable("Bending Enabled", params_.bendingEnabled);

    viewer.ngui->addWindow(Eigen::Vector2i(1000, 0), "New Objects");
    viewer.ngui->addGroup("New Particles");
    viewer.ngui->addVariable("Is Fixed", params_.particleFixed);
    viewer.ngui->addVariable("Mass", params_.particleMass);
    viewer.ngui->addGroup("New Saws");
    viewer.ngui->addVariable("Radius", params_.sawRadius);
    viewer.ngui->addGroup("New Springs");
    viewer.ngui->addVariable("Max Spring Dist", params_.maxSpringDist);
    viewer.ngui->addVariable("Base Stiffness", params_.springStiffness);
    viewer.ngui->addGroup("New Rods");
    viewer.ngui->addVariable("Num Segments", params_.rodSegments);
    viewer.ngui->addVariable("Density", params_.rodDensity);
    viewer.ngui->addVariable("Stretching Stiffness", params_.rodStretchingStiffness);
    viewer.ngui->addVariable("Bending Stiffness", params_.rodBendingStiffness);

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

    
    for(std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
    {
        switch((*it)->getType())
        {
        case SimParameters::CT_SPRING:
        {
            Eigen::Vector3d color;
            if((*it)->associatedBendingStencils.empty())
                color << 0.0, 0.0, 1.0;
            else
                color << 0.75, 0.5, 0.75;
            Vector2d sourcepos = particles_[(*it)->p1].pos;
            Vector2d destpos   = particles_[(*it)->p2].pos;

            Vector2d vec = destpos - sourcepos;
            Vector2d perp(-vec[1], vec[0]);
            perp /= perp.norm();

            double dist = (sourcepos-destpos).norm();

            double width = baselinewidth/(1.0+ 20.0 * dist * dist);

            for (int i = 0; i < 4; i++)
                vertexColors.push_back(color);

            verts.push_back(Eigen::Vector3d(sourcepos[0] + width * perp[0], sourcepos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(sourcepos[0] - width * perp[0], sourcepos[1] - width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] + width * perp[0], destpos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] - width * perp[0], destpos[1] - width * perp[1], -eps));

            faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));
            faces.push_back(Eigen::Vector3i(idx + 2, idx + 1, idx + 3));
            idx += 4;
            
            break;
        }
        case SimParameters::CT_RIGIDROD:
        {
            Eigen::Vector3d color;
            if((*it)->associatedBendingStencils.empty())
                color << 1.0, 0.0, 1.0;
            else
                color << 1.0, 1.0, 0.3;

            Vector2d sourcepos = particles_[(*it)->p1].pos;
            Vector2d destpos   = particles_[(*it)->p2].pos;
            Vector2d vec = destpos - sourcepos;
            Vector2d perp(-vec[1], vec[0]);
            perp /= perp.norm();

            double width = baselinewidth;

            for (int i = 0; i < 4; i++)
                vertexColors.push_back(color);

            verts.push_back(Eigen::Vector3d(sourcepos[0] + width * perp[0], sourcepos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(sourcepos[0] - width * perp[0], sourcepos[1] - width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] + width * perp[0], destpos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] - width * perp[0], destpos[1] - width * perp[1], -eps));

            faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));
            faces.push_back(Eigen::Vector3i(idx + 2, idx + 1, idx + 3));
            idx += 4;

            break;
        }
        default:
            break;
        }
    }

    int nparticles = particles_.size();

    for(int i=0; i<nparticles; i++)
    {
        double radius = baseradius*sqrt(getTotalParticleMass(i));
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
    bendingStencils_.clear();
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
    // TODO handle constraints and flexible rods
    VectorXd q, v;
    buildConfiguration(q, v);
    numericalIntegration(q, v);
    unbuildConfiguration(q, v);

    pruneOverstrainedSprings();
    deleteSawedObjects();
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

    int numparticles = particles_.size()-1;

    for(int i=0; i<numparticles; i++)
    {
        if(particles_[i].inert)
            continue;
        Vector2d pos = particles_[i].pos;
        double dist = (pos-newpos).norm();
        if(dist <= params_.maxSpringDist)
        {
            switch(params_.connectorType)
            {
            case SimParameters::CT_SPRING:
            {
                connectors_.push_back(new Spring(newid, i, 0, params_.springStiffness/dist, dist, true));
                break;
            }
            case SimParameters::CT_RIGIDROD:
            {
                connectors_.push_back(new RigidRod(newid, i, 0, dist));
                break;
            }
            case SimParameters::CT_FLEXROD:
            {
                // TODO add flexible rods here
                break;
            }
            default:
                break;
            }
        }
    }
}

double GooHook::getTotalParticleMass(int idx)
{
    double mass = particles_[idx].mass;
    //TODO account for flexible rods adding mass
    return mass;
}

void GooHook::addSaw(double x, double y)
{
    saws_.push_back(Saw(Vector2d(x,y), params_.sawRadius));    
}

void GooHook::buildConfiguration(VectorXd &q, VectorXd &v)
{
    int ndofs = 2*particles_.size();
    q.resize(ndofs);
    v.resize(ndofs);

    for(int i=0; i<(int)particles_.size(); i++)
    {
        q.segment<2>(2*i) = particles_[i].pos;
        v.segment<2>(2*i) = particles_[i].vel;
    }
}

void GooHook::unbuildConfiguration(const VectorXd &q, const VectorXd &v)
{
    int ndofs = q.size();
    assert(ndofs == int(2*particles_.size()));

    for(int i=0; i<ndofs/2; i++)
    {
        particles_[i].pos = q.segment<2>(2*i);
        particles_[i].vel = v.segment<2>(2*i);
    }
}

int GooHook::getNumRigidRods()
{
    int nrods = 0;
    for(std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
    {
        if((*it)->getType() == SimParameters::CT_RIGIDROD)
            nrods++;
    }
    return nrods;
}

void GooHook::numericalIntegration(VectorXd &q, VectorXd &v)
{
    // TODO handle constraints and flexible rods
    VectorXd F;
    SparseMatrix<double> H;
    SparseMatrix<double> Minv;

    computeMassInverse(Minv);

    VectorXd oldq = q;

    q += params_.timeStep*v;
    computeForceAndHessian(q, oldq, F, H);

    if (params_.constraintHandling == SimParameters::CH_PENALTY) {
        v += params_.timeStep*Minv*F;
    } else if (params_.constraintHandling == SimParameters::CH_STEPPROJECT) {
        v += params_.timeStep*Minv*F;
        std::vector<RigidRod *> rods;
        int nsprings = (int)connectors_.size();
        int nparticles = (int)particles_.size();
        int m_start = 2 * nparticles;

        for(int i=0; i<nsprings; i++)
        {
            if(connectors_[i]->getType() != SimParameters::CT_RIGIDROD)
                continue;

            RigidRod *r = (RigidRod *)connectors_[i];
            rods.push_back(r);
        }

        int x_size = m_start + int(rods.size());
        VectorXd f(x_size);
        std::vector<Eigen::Triplet<double> > df_coeffs;

        VectorXd q0 = q;
        processRodProjection(rods, q, q0, f, df_coeffs);

        int nIters = 0;
        while (f.norm() > params_.NewtonTolerance) {
            SparseMatrix<double> df(x_size, x_size);
            df.setFromTriplets(df_coeffs.begin(), df_coeffs.end());
            Eigen::SparseQR<Eigen::SparseMatrix<double>, COLAMDOrdering<int>> solver;
            solver.compute(df);

            VectorXd delta = solver.solve(-f);
            q0 += delta.segment(0, 2*nparticles);
            for (int i = 0; i < rods.size(); i++) {
                RigidRod *r = rods[i];
                r->lambda += delta(m_start + i);
            }
            processRodProjection(rods, q, q0, f, df_coeffs);

            nIters++;
            std::cout << "+++ " << nIters << " iter +++\n\n\n";
            if (nIters >= params_.NewtonMaxIters) {
                break;
            }
        }

        v += (q0 - q) / params_.timeStep;
        q = q0;
    } else {
        std::vector<RigidRod *> rods;
        int nsprings = (int)connectors_.size();
        int nparticles = (int)particles_.size();

        for(int i=0; i<nsprings; i++)
        {
            if(connectors_[i]->getType() != SimParameters::CT_RIGIDROD)
                continue;

            RigidRod *r = (RigidRod *)connectors_[i];
            rods.push_back(r);
        }

        int m = rods.size();
        VectorXd lambda(m);
        for (int i = 0; i < m; i++) {
            lambda[i] = rods[i]->lambda;
        }

        VectorXd q2_cons = q + params_.timeStep*v +
                             params_.timeStep*params_.timeStep*Minv*F;

        SparseMatrix<double> dg(m, 2 * nparticles);
        std::vector<Eigen::Triplet<double> > dg_coeffs;
        compute_dg(rods, q, dg_coeffs);
        dg.setFromTriplets(dg_coeffs.begin(), dg_coeffs.end());
        // std::cout << "dg=[\n" << dg << "]" << std::endl;
        SparseMatrix<double> dgT = dg.transpose();
        SparseMatrix<double> dq2 = params_.timeStep*params_.timeStep*Minv*dgT;

        VectorXd q2_0 = q2_cons + dq2*lambda;
        // std::cout << "q_i+2=[\n" << q2_0 << "]" << std::endl;
        VectorXd f(m);
        compute_g(rods, q2_0, f);
        // std::cout << "f=[\n" << f << "]" << std::endl;

        int nIters = 0;
        while (f.norm() > params_.NewtonTolerance) {
            SparseMatrix<double> dg_q2(m, 2 * nparticles);
            std::vector<Eigen::Triplet<double> > dg_q2_coeffs;
            compute_dg(rods, q2_0, dg_q2_coeffs);
            dg_q2.setFromTriplets(dg_q2_coeffs.begin(), dg_q2_coeffs.end());
            // std::cout << "dg(q_i+2)=[\n" << dg_q2 << "]" << std::endl;
            SparseMatrix<double> df = dg_q2 * dq2;

            Eigen::SparseQR<Eigen::SparseMatrix<double>, COLAMDOrdering<int>> solver;
            solver.compute(df);
            VectorXd delta = solver.solve(-f);
            lambda += delta;
            q2_0 = q2_cons + dq2*lambda;
            compute_g(rods, q2_0, f);

            nIters++;
            // std::cout << "+++ " << nIters << " iter +++\n\n\n";
            if (nIters >= params_.NewtonMaxIters) {
                break;
            }
        }

        v += params_.timeStep*Minv*F + params_.timeStep*Minv*dgT*lambda;
        for (int i = 0; i < m; i++) {
            rods[i]->lambda = lambda[i];
        }
    }
}

void GooHook::computeMassInverse(Eigen::SparseMatrix<double> &Minv)
{
    int ndofs = 2*int(particles_.size());

    Minv.resize(ndofs, ndofs);
    Minv.setZero();

    std::vector<Eigen::Triplet<double> > Minvcoeffs;
    for(int i=0; i<ndofs/2; i++)
    {
        Minvcoeffs.push_back(Eigen::Triplet<double>(2*i,   2*i,   1.0/getTotalParticleMass(i)));
        Minvcoeffs.push_back(Eigen::Triplet<double>(2*i+1, 2*i+1, 1.0/getTotalParticleMass(i)));
    }

    Minv.setFromTriplets(Minvcoeffs.begin(), Minvcoeffs.end());
}

void GooHook::computeForceAndHessian(const VectorXd &q, const VectorXd &qprev, Eigen::VectorXd &F, SparseMatrix<double> &H)
{
    F.resize(q.size());
    F.setZero();
    H.resize(q.size(), q.size());
    H.setZero();

    std::vector<Eigen::Triplet<double> > Hcoeffs;
    if(params_.gravityEnabled)
        processGravityForce(F);
    if(params_.springsEnabled)
        processSpringForce(q, F, Hcoeffs);
    if(params_.dampingEnabled)
        processDampingForce(q, qprev, F, Hcoeffs);
    if(params_.floorEnabled)
        processFloorForce(q, qprev, F, Hcoeffs);
    if(params_.constraintHandling == SimParameters::CH_PENALTY)
        processRodPenaltyForce(q, F);

    H.setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());
}

void GooHook::processGravityForce(VectorXd &F)
{
    int nparticles = (int)particles_.size();
    for(int i=0; i<nparticles; i++)
    {
        if(!particles_[i].fixed)
        {
            F[2*i+1] += params_.gravityG*getTotalParticleMass(i);
        }
    }
}

void GooHook::processRodPenaltyForce(const Eigen::VectorXd &q, Eigen::VectorXd &F)
{
    int nsprings = (int)connectors_.size();

    for(int i=0; i<nsprings; i++)
    {
        if(connectors_[i]->getType() != SimParameters::CT_RIGIDROD)
            continue;

        RigidRod &r = *(RigidRod *)connectors_[i];
        Vector2d p1 = q.segment<2>(2*r.p1);
        Vector2d p2 = q.segment<2>(2*r.p2);
        Vector2d diff = p1 - p2;
        double g = diff.squaredNorm() - r.length * r.length;
        Vector2d localF = -4 * params_.penaltyStiffness * g * diff;
        F.segment<2>(2*r.p1) += localF;
        F.segment<2>(2*r.p2) -= localF;
    }
}

// dg(q) is mx2n
void GooHook::compute_dg(const std::vector<RigidRod *> &rods,
                                 const Eigen::VectorXd &q, std::vector<Eigen::Triplet<double> > &dg)
{
    int nrods = (int)rods.size();
    dg.clear();
    for(int i=0; i<nrods; i++)
    {
        RigidRod &r = *rods[i];
        Vector2d p1 = q.segment<2>(2*r.p1);
        Vector2d p2 = q.segment<2>(2*r.p2);
        Vector2d diff = p1 - p2;
        Vector2d local_dg = 2 * diff;

        dg.push_back(Eigen::Triplet<double>(i, 2*r.p1, local_dg[0]));
        dg.push_back(Eigen::Triplet<double>(i, 2*r.p1+1, local_dg[1]));
        dg.push_back(Eigen::Triplet<double>(i, 2*r.p2, -local_dg[0]));
        dg.push_back(Eigen::Triplet<double>(i, 2*r.p2+1, -local_dg[1]));
    }
}

void GooHook::compute_g(const std::vector<RigidRod *> &rods,
                                 const Eigen::VectorXd &q, Eigen::VectorXd &g)
{
    int nrods = (int)rods.size();
    for(int i=0; i<nrods; i++)
    {
        RigidRod &r = *rods[i];
        Vector2d p1 = q.segment<2>(2*r.p1);
        Vector2d p2 = q.segment<2>(2*r.p2);
        Vector2d diff = p1 - p2;
        double local_g = diff.squaredNorm() - r.length * r.length;

        g[i] = local_g;
    }
}

void GooHook::processRodProjection(const std::vector<RigidRod *> &rods,
                                   const Eigen::VectorXd &q_telda, const Eigen::VectorXd &q,
                                   Eigen::VectorXd &f, std::vector<Eigen::Triplet<double> > &df)
{
    int nparticles = (int)particles_.size();
    int nrods = (int)rods.size();
    int m = nrods;
    SparseMatrix<double> Minv;
    computeMassInverse(Minv);

    f.setZero();
    VectorXd sigma(2 * nparticles);
    sigma.setZero();
    int m_start = 2 * nparticles;
    for(int i=0; i<nrods; i++)
    {
        RigidRod &r = *rods[i];
        Vector2d p1 = q.segment<2>(2*r.p1);
        Vector2d p2 = q.segment<2>(2*r.p2);
        Vector2d diff = p1 - p2;
        double g = diff.squaredNorm() - r.length * r.length;

        VectorXd dgi(2 * nparticles);
        dgi.setZero();
        dgi.segment<2>(2*r.p1) = 2 * diff;
        dgi.segment<2>(2*r.p2) = -2 * diff;
        sigma += r.lambda * Minv * dgi;

        f(m_start + i) = g;
    }
    f.segment(0, 2 * nparticles) = q - q_telda + sigma;

    df.clear();
    for (int i = 0; i < nparticles; i++) {
        df.push_back(Eigen::Triplet<double>(2*i, 2*i, 1.0));
        df.push_back(Eigen::Triplet<double>(2*i+1, 2*i+1, 1.0));
    }

    for(int i=0; i<nrods; i++)
    {
        RigidRod &r = *rods[i];
        Vector2d p1 = q.segment<2>(2*r.p1);
        Vector2d p2 = q.segment<2>(2*r.p2);
        Vector2d diff = p1 - p2;
        double l = r.lambda;

        // Push M_inv * dg_{i}^T
        Vector2d local_dg = 2 * diff;
        double m1 = Minv.coeff(2*r.p1, 2*r.p1);
        double m2 = Minv.coeff(2*r.p2, 2*r.p2);
        df.push_back(Eigen::Triplet<double>(2*r.p1, m_start + i, m1 * local_dg[0]));
        df.push_back(Eigen::Triplet<double>(2*r.p1 + 1, m_start + i, m1 * local_dg[1]));
        df.push_back(Eigen::Triplet<double>(2*r.p2, m_start + i, -m2 * local_dg[0]));
        df.push_back(Eigen::Triplet<double>(2*r.p2 + 1, m_start + i, -m2 * local_dg[1]));

        // Push dgi
        df.push_back(Eigen::Triplet<double>(m_start + i, 2*r.p1, local_dg[0]));
        df.push_back(Eigen::Triplet<double>(m_start + i, 2*r.p1+1, local_dg[1]));
        df.push_back(Eigen::Triplet<double>(m_start + i, 2*r.p2, -local_dg[0]));
        df.push_back(Eigen::Triplet<double>(m_start + i, 2*r.p2+1, -local_dg[1]));

        // Push lambda * M_inv * Heissan
        df.push_back(Eigen::Triplet<double>(2*r.p1, 2*r.p1, 2 * m1 * l));
        df.push_back(Eigen::Triplet<double>(2*r.p1+1, 2*r.p1+1, 2 * m1 * l));
        df.push_back(Eigen::Triplet<double>(2*r.p1, 2*r.p2, -2 * m1 * l));
        df.push_back(Eigen::Triplet<double>(2*r.p1+1, 2*r.p2+1, -2 * m1 * l));
        df.push_back(Eigen::Triplet<double>(2*r.p2, 2*r.p1, -2 * m1 * l));
        df.push_back(Eigen::Triplet<double>(2*r.p2+1, 2*r.p1+1, -2 * m1 * l));
        df.push_back(Eigen::Triplet<double>(2*r.p2, 2*r.p2, 2 * m1 * l));
        df.push_back(Eigen::Triplet<double>(2*r.p2+1, 2*r.p2+1, 2 * m1 * l));
    }
}

void GooHook::processSpringForce(const Eigen::VectorXd &q, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H)
{
    int nsprings = (int)connectors_.size();

    for(int i=0; i<nsprings; i++)
    {
        if(connectors_[i]->getType() != SimParameters::CT_SPRING)
            continue;
        Spring &s = *(Spring *)connectors_[i];
        Vector2d p1 = q.segment<2>(2*s.p1);
        Vector2d p2 = q.segment<2>(2*s.p2);
        double dist = (p2-p1).norm();
        Vector2d localF = s.stiffness*(dist-s.restlen)/dist * (p2-p1);
        F.segment<2>(2*s.p1) += localF;
        F.segment<2>(2*s.p2) -= localF;

        Matrix2d I;
        I << 1, 0, 0, 1;
        Matrix2d localH = s.stiffness * (1.0 - s.restlen/dist)*I;
        localH += s.stiffness*s.restlen*(p2-p1)*(p2-p1).transpose()/dist/dist/dist;

        for(int j=0; j<2; j++)
            for(int k=0; k<2;k++)
            {
                H.push_back(Eigen::Triplet<double>(2*s.p1+j, 2*s.p1+k, localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p2+j, 2*s.p2+k, localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p1+j, 2*s.p2+k, -localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p2+j, 2*s.p1+k, -localH.coeff(j,k)));
            }
    }
}

void GooHook::processDampingForce(const VectorXd &q, const VectorXd &qprev, VectorXd &F, std::vector<Eigen::Triplet<double> > &H)
{
    int nsprings = (int)connectors_.size();

    for(int i=0; i<nsprings; i++)
    {
        if(connectors_[i]->getType() != SimParameters::CT_SPRING)
            continue;
        Spring &s = *(Spring *)connectors_[i];
        Vector2d p1 = q.segment<2>(2*s.p1);
        Vector2d p2 = q.segment<2>(2*s.p2);
        Vector2d p1prev = qprev.segment<2>(2*s.p1);
        Vector2d p2prev = qprev.segment<2>(2*s.p2);

        Vector2d relvel = (p2 - p2prev)/params_.timeStep - (p1 - p1prev)/params_.timeStep;
        Vector2d localF = params_.dampingStiffness*relvel;
        F.segment<2>(2*s.p1) += localF;
        F.segment<2>(2*s.p2) -= localF;

        Matrix2d I;
        I << 1, 0, 0, 1;
        Matrix2d localH = params_.dampingStiffness*I/params_.timeStep;

        for(int j=0; j<2; j++)
            for(int k=0; k<2;k++)
            {
                H.push_back(Eigen::Triplet<double>(2*s.p1+j, 2*s.p1+k, localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p2+j, 2*s.p2+k, localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p1+j, 2*s.p2+k, -localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p2+j, 2*s.p1+k, -localH.coeff(j,k)));
            }
    }
}

void GooHook::processFloorForce(const VectorXd &q, const VectorXd &qprev, VectorXd &F, std::vector<Eigen::Triplet<double> > &H)
{
    int nparticles = particles_.size();

    double basestiffness = 10000;
    double basedrag = 1000.0;

    for(int i=0; i<nparticles; i++)
    {
        if(q[2*i+1] < -0.5 && ! particles_[i].fixed)
        {
            double vel = (q[2*i+1]-qprev[2*i+1])/params_.timeStep;
            double dist = -0.5 - q[2*i+1];

            F[2*i+1] += basestiffness*dist - basedrag*dist*vel;

            H.push_back(Eigen::Triplet<double>(2*i+1, 2*i+1, basestiffness
                - 0.5*basedrag/params_.timeStep
                + basedrag*qprev[2*i+1]/params_.timeStep
                - 2.0*basedrag*q[2*i+1]/params_.timeStep));
        }
    }
}

double GooHook::ptSegmentDist(const Vector2d &p, const Vector2d &q1, const Vector2d &q2)
{
    double t = (p-q1).dot(q2-q1) / (q2-q1).dot(q2-q1);
    double linedistsq = (q1 + t*(q2-q1) - p).squaredNorm();
    double q1dist = (p-q1).squaredNorm();
    double q2dist = (p-q2).squaredNorm();
    double mindistsq = std::min(linedistsq, std::min(q1dist, q2dist));
    return sqrt(mindistsq);
}

void GooHook::detectSawedConnectors(std::set<int> &connectorsToDelete)
{
    for(int i=0; i<(int)connectors_.size(); i++)
    {
        Vector2d pos1 = particles_[connectors_[i]->p1].pos;
        Vector2d pos2 = particles_[connectors_[i]->p2].pos;
        double maxx = std::max(pos1[0], pos2[0]);
        double minx = std::min(pos1[0], pos2[0]);
        double maxy = std::max(pos1[1], pos2[1]);
        double miny = std::min(pos1[1], pos2[1]);
        for(std::vector<Saw>::iterator saw = saws_.begin(); saw != saws_.end(); ++saw)
        {
            Vector2d sawpos = saw->pos;
            double sawr = saw->radius;

            if(sawpos[0] - sawr > maxx || sawpos[0] + sawr < minx || sawpos[1] - sawr > maxy || sawpos[1] + sawr < miny)
                continue;

            double sawspringdist = ptSegmentDist(sawpos, pos1, pos2);
            if(sawspringdist <= sawr)
            {
                connectorsToDelete.insert(i);
                break;
            }
        }
    }
}

void GooHook::detectSawedParticles(std::set<int> &particlesToDelete)
{
    for(int i=0; i<(int)particles_.size(); i++)
    {
        Vector2d partpos = particles_[i].pos;

        if(!(fabs(partpos[0]) < 2 && fabs(partpos[1]) < 2))
        {
            particlesToDelete.insert(i);
            break;
        }

        for(std::vector<Saw>::iterator it = saws_.begin(); it != saws_.end(); ++it)
        {
            Vector2d sawpos = it->pos;
            double sqdist = (sawpos-partpos).squaredNorm();
            if(sqdist < it->radius*it->radius)
            {
                particlesToDelete.insert(i);
                break;
            }
        }
    }
}

void GooHook::deleteSawedObjects()
{
    std::set<int> particlestodelete;
    std::set<int> connectorstodelete;
    std::set<int> bendingtodelete;
    detectSawedParticles(particlestodelete);
    detectSawedConnectors(connectorstodelete);

    std::vector<Particle, Eigen::aligned_allocator<Particle>> newparticles;
    std::vector<Connector *> newconnectors;
    std::vector<BendingStencil> newbending;
    std::vector<int> remainingparticlemap;
    std::vector<int> remainingbendingmap;

    if(!particlestodelete.empty())
    {
        for(int i=0; i<(int)connectors_.size(); i++)
        {
            if(particlestodelete.count(connectors_[i]->p1) || particlestodelete.count(connectors_[i]->p2))
                connectorstodelete.insert(i);
        }

        for(int i=0; i<(int)particles_.size(); i++)
        {
            if(particlestodelete.count(i) == 0)
            {
                remainingparticlemap.push_back(newparticles.size());
                newparticles.push_back(particles_[i]);
            }
            else
                remainingparticlemap.push_back(-1);
        }
    }
    if(!connectorstodelete.empty())
    {
        for(std::set<int>::iterator it = connectorstodelete.begin(); it != connectorstodelete.end(); ++it)
        {
            for(std::set<int>::iterator bit = connectors_[*it]->associatedBendingStencils.begin(); bit != connectors_[*it]->associatedBendingStencils.end(); ++bit)
            {
                bendingtodelete.insert(*bit);
            }
        }
        for(int i=0; i<(int)connectors_.size(); i++)
        {
            if(connectorstodelete.count(i) == 0)
            {
                newconnectors.push_back(connectors_[i]);
            }
            else
                delete connectors_[i];
        }
    }
    if(!bendingtodelete.empty())
    {
        int newidx=0;
        for(int i=0; i<(int)bendingStencils_.size(); i++)
        {
            if(bendingtodelete.count(i) == 0)
            {
                newbending.push_back(bendingStencils_[i]);
                remainingbendingmap.push_back(newidx++);
            }
            else
                remainingbendingmap.push_back(-1);
        }
    }

    if (!connectorstodelete.empty() || !particlestodelete.empty())
    {
        if (!connectorstodelete.empty())
            connectors_ = newconnectors;
        if (!bendingtodelete.empty())
        {
            bendingStencils_ = newbending;
            for (std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
            {
                std::set<int> newass;
                for (std::set<int>::iterator sit = (*it)->associatedBendingStencils.begin(); sit != (*it)->associatedBendingStencils.end(); ++sit)
                {
                    if (bendingtodelete.count(*sit) == 0)
                        newass.insert(remainingbendingmap[*sit]);
                }
                (*it)->associatedBendingStencils = newass;
            }
        }
        if (!particlestodelete.empty())
        {
            particles_ = newparticles;
            for (std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
            {
                (*it)->p1 = remainingparticlemap[(*it)->p1];
                (*it)->p2 = remainingparticlemap[(*it)->p2];
            }
            for (std::vector<BendingStencil>::iterator it = bendingStencils_.begin(); it != bendingStencils_.end(); ++it)
            {
                it->p1 = remainingparticlemap[it->p1];
                it->p2 = remainingparticlemap[it->p2];
                it->p3 = remainingparticlemap[it->p3];
            }
        }
    }
}

void GooHook::pruneOverstrainedSprings()
{
    int nsprings = connectors_.size();

    std::vector<int> toremove;
    for (int i = 0; i < nsprings; i++)
    {
        if (connectors_[i]->getType() != SimParameters::CT_SPRING)
            continue;

        Spring &s = *(Spring *)connectors_[i];
        if (s.canSnap)
        {
            Vector2d srcpos = particles_[s.p1].pos;
            Vector2d dstpos = particles_[s.p2].pos;
            double dist = (dstpos - srcpos).norm();

            double strain = (dist - s.restlen) / s.restlen;
            if (strain > params_.maxSpringStrain)
                toremove.push_back(i);
        }
    }

    for (std::vector<int>::reverse_iterator it = toremove.rbegin(); it != toremove.rend(); ++it)
    {
        assert(connectors_[*it]->associatedBendingStencils.empty());
        delete connectors_[*it];
        connectors_.erase(connectors_.begin() + *it);
    }
}
