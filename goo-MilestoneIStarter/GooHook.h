#include "PhysicsHook.h"
#include "SceneObjects.h"
#include <deque>
#include <map>
#include "SimParameters.h"
#include <Eigen/Sparse>
#include <Eigen/StdVector>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

struct MouseClick
{
    double x;
    double y;
    SimParameters::ClickMode mode;
};

class GooHook : public PhysicsHook
{
public:
    GooHook() : PhysicsHook() {}

    virtual void initGUI(igl::viewer::Viewer &viewer);

    virtual void initSimulation();

    virtual void mouseClicked(double x, double y, int button)
    {
        message_mutex.lock();
        {
            MouseClick mc;
            mc.x = x;
            mc.y = y;
            mc.mode = params_.clickMode;
            mouseClicks_.push_back(mc);
        }
        message_mutex.unlock();
    }

    virtual void updateRenderGeometry();
    
    virtual void tick();

    virtual bool simulateOneStep();

    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer)
    {
        viewer.data.clear();
        viewer.data.set_mesh(renderQ, renderF);
        viewer.data.set_colors(renderC);
    }

private:
    SimParameters params_;
    double time_;
    std::vector<Particle, Eigen::aligned_allocator<Particle> > particles_;
    std::vector<Connector *> connectors_;
    std::vector<Saw> saws_;
    
    std::mutex message_mutex;
    std::deque<MouseClick> mouseClicks_;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd renderC;

    void addParticle(double x, double y);
    void addSaw(double x, double y);      

    void removeObjects();

    Eigen::VectorXd configVector(); 
    Eigen::VectorXd prevConfigVector(); 
    Eigen::VectorXd configVelVector();
    Eigen::VectorXd gravity();
    Eigen::MatrixXd gravityHeissan();
    Eigen::VectorXd springForce(Eigen::VectorXd q);
    Eigen::VectorXd springForceHeissan(Eigen::VectorXd q);
    Eigen::VectorXd viscousDamping(Eigen::VectorXd v);
    Eigen::VectorXd viscousDampingHeissan(Eigen::VectorXd v);
    // V = Kmg(y + 0.5) if y < -0.5
    // V = 0 if y >= -0.5
    Eigen::VectorXd floorForce(Eigen::VectorXd q, Eigen::VectorXd v);
    Eigen::MatrixXd floorForceHeissan(Eigen::VectorXd q, Eigen::VectorXd v);
    SpMat selector(int i);

    SpMat massMatrix();
    SpMat massInvMatrix();
};
