#pragma once

#include <algorithm>
#include <math.h>
#include <stdio.h>

#include "Solid.h"

using namespace std;
using namespace Eigen;

class MACGrid {
    VectorXd *src_;
    VectorXd *dst_;
    VectorXd normalX_;
    VectorXd normalY_;
    VectorXi cell_;
    VectorXi body_;
    VectorXi mask_;

    int w_;
    int h_;

    /* X and Y offset from top left grid cell.
     * This is (0.5,0.5) for centered quantities such as density,
     * and (0.0, 0.5) or (0.5, 0.0) for jittered quantities like the velocity.
     */
    double ox_;
    double oy_;

    double hx_;

    double lerp(double a, double b, double x) const {
        return a * (1.0 - x) + b * x;
    }
    
    Vector2d back_advect(double x, double y, double dt, const MACGrid &u, const MACGrid &v) const {
        // std::cout << "back_advect: " << x << "," << y << "\n";
        double vx = u.lerp(x, y) / hx_;
        double vy = v.lerp(x, y) / hx_;
        
        x -= vx * dt;
        y -= vy * dt;
        return Vector2d(x, y);
    }
    
public:
    MACGrid(int w, int h, double ox, double oy, double hx)
            : w_(w), h_(h), ox_(ox), oy_(oy), hx_(hx) {
        src_ = new VectorXd(w_ * h_);
        dst_ = new VectorXd(w_ * h_);
                
        src_->setZero();
        dst_->setZero();

        normalX_.resize(w_ * h_);
        normalY_.resize(w_ * h_);

        cell_.resize(w_ * h_);
        body_.resize(w_ * h_);
        mask_.resize(w_ * h_);

        cell_.setZero();
    }
    
    ~MACGrid() {
        delete src_;
        delete dst_;
    }

    const VectorXi &cell() const {
        return cell_;
    }

    const VectorXi &body() const {
        return body_;
    }
    
    void flip() {
        swap(src_, dst_);
    }
    
    const VectorXd &src() const {
        return *src_;
    }
    
    double at(int x, int y) const {
        return (*src_)(x + y * w_);
    }
    
    double &at(int x, int y) {
        return (*src_)(x + y * w_);
    }
    
    double lerp(double x, double y) const {
        x = min(max(x - ox_, 0.0), w_ - 1.001);
        y = min(max(y - oy_, 0.0), h_ - 1.001);
        int ix = (int)x;
        int iy = (int)y;
        double dx = x - ix;
        double dy = y - iy;
        
        double x00 = at(ix + 0, iy + 0), x10 = at(ix + 1, iy + 0);
        double x01 = at(ix + 0, iy + 1), x11 = at(ix + 1, iy + 1);
        // std::cout << "lerp: " << ix << "," << iy << "\n";
        // std::cout << "lerp: " << x00 << "," << x01 << "," << x10 << "," << x11 << "\n";
        return lerp(lerp(x00, x10, dx), lerp(x01, x11, dx), dy);
    }


    
    /* If the point (x, y) is inside a solid, project it back out to the
     * closest point on the surface of the solid.
     */
    Vector2d backProject(double x, double y, const vector<SolidBody *> &bodies) {
        int rx = min(max((int)(x - ox_), 0), w_ - 1);
        int ry = min(max((int)(y - oy_), 0), h_ - 1);
        
        if (cell_(rx + ry * w_) != CELL_FLUID) {
            x = (x - ox_) * hx_;
            y = (y - oy_) * hx_;
            bodies[body_(rx + ry * w_)]->closestSurfacePoint(x, y);
            x = x / hx_ + ox_;
            y = y / hx_ + oy_;
        }

        return Vector2d(x, y);
    }

    void advect(double dt, const MACGrid &u, const MACGrid &v, const vector<SolidBody *> &bodies) {
        for (int iy = 0; iy < h_; iy++) {
            for (int ix = 0; ix < w_; ix++) {
                double x = ix + ox_;
                double y = iy + oy_;
                
                auto xo = back_advect(x, y, dt, u, v);
                // std::cout << "xo = [" << xo << "]\n";
                xo = backProject(xo(0), xo(1), bodies);
                (*dst_)(ix + iy * w_) = lerp(xo(0), xo(1));
            }
        }
    }

    void addInflow(double x0, double y0, double x1, double y1, double v) {
        int ix0 = (int)(x0/hx_ - ox_);
        int iy0 = (int)(y0/hx_ - oy_);
        int ix1 = (int)(x1/hx_ - ox_);
        int iy1 = (int)(y1/hx_ - oy_);
        
        for (int y = max(iy0, 0); y < min(iy1, h_); y++) {
            for (int x = max(ix0, 0); x < min(ix1, h_); x++) {
                (*src_)(x + y * w_) = v;
            }
        }
    }

    
    /* Fill all solid related fields - that is, _cell, _body and _normalX/Y */
    void fillSolidFields(const vector<SolidBody *> &bodies) {
        if (bodies.empty())
            return;
        
        for (int iy = 0, idx = 0; iy < h_; iy++) {
            for (int ix = 0; ix < w_; ix++, idx++) {
                double x = (ix + ox_) * hx_;
                double y = (iy + oy_) * hx_;
                
                /* Search closest solid */
                body_(idx) = 0;
                double d = bodies[0]->distance(x, y);
                for (unsigned i = 1; i < bodies.size(); i++) {
                    double id = bodies[i]->distance(x, y);
                    if (id < d) {
                        body_(idx) = i;
                        d = id;
                    }
                }
                
                /* If distance to closest solid is negative, this cell must be
                 * inside it
                 */
                if (d < 0.0)
                    cell_(idx) = CELL_SOLID;
                else
                    cell_(idx) = CELL_FLUID;
                double nx, ny;
                bodies[body_(idx)]->distanceNormal(nx, ny, x, y);
                normalX_(idx) = nx;
                normalY_(idx) = ny;
            }
        }
    }

    void fillSolidMask() {
        for (int y = 1; y < h_ - 1; y++) {
            for (int x = 1; x < w_ - 1; x++) {
                int idx = x + y * w_;
                
                if (cell_(idx) == CELL_FLUID)
                    continue;
                
                double nx = normalX_(idx);
                double ny = normalY_(idx);

                mask_(idx) = 0;
                if (nx != 0.0 && cell_(idx + sgn(nx))    != CELL_FLUID)
                    mask_(idx) |= 1; /* Neighbour in normal x direction is blocked */
                if (ny != 0.0 && cell_(idx + sgn(ny) * w_) != CELL_FLUID)
                    mask_(idx) |= 2; /* Neighbour in normal y direction is blocked */
            }
        }
    }

   /* Solve for value at index idx using values of neighbours in normal x/y
     * direction. The value is computed such that the directional derivative
     * along distance field normal is 0.
     */
    double extrapolateNormal(int idx) {
        double nx = normalX_(idx);
        double ny = normalY_(idx);
        
        double srcX = (*src_)(idx + sgn(nx));
        double srcY = (*src_)(idx + sgn(ny) * w_);
        
        return (fabs(nx)*srcX + fabs(ny)*srcY)/(fabs(nx) + fabs(ny));
    }
    
    /* Given that a neighbour in upstream direction specified by mask (1=x, 2=y)
     * now has been solved for, update the mask appropriately and, if this cell
     * can now be computed, add it to the queue of ready cells
     */
    void freeNeighbour(int idx, stack<int> &border, int mask) {
        mask_(idx) &= ~mask;
        if (cell_(idx) != CELL_FLUID && mask_(idx) == 0)
            border.push(idx);
    }
    
    void extrapolate() {
        fillSolidMask();

        /* Queue of cells which can be computed */
        stack<int> border;
        /* Initialize queue by finding all solid cells with mask=0 (ready for
         * extrapolation)
         */
        for (int y = 1; y < h_ - 1; y++) {
            for (int x = 1; x < w_ - 1; x++) {
                int idx = x + y * w_;

                if (cell_(idx) != CELL_FLUID && mask_(idx) == 0)
                    border.push(idx);
            }
        }

        while (!border.empty()) {
            int idx = border.top();
            border.pop();

            /* Solve for value in cell */
            (*src_)(idx) = extrapolateNormal(idx);

            /* Notify adjacent cells that this cell has been computed and can
             * be used as an upstream value
             */
            if (normalX_(idx - 1) > 0.0)
                freeNeighbour(idx -  1, border, 1);
            if (normalX_(idx + 1) < 0.0)
                freeNeighbour(idx +  1, border, 1);
            if (normalY_(idx - w_) > 0.0)
                freeNeighbour(idx - w_, border, 2);
            if (normalY_(idx + w_) < 0.0)
                freeNeighbour(idx + w_, border, 2);
        }
    }
};

struct BodyForces {
    double g;
};

/* Fluid solver class. Sets up the fluid quantities, forces incompressibility
 * performs advection and adds inflows.
 */
class FluidSolver {
    /* Fluid quantities */
    MACGrid *d_;
    MACGrid *u_;
    MACGrid *v_;
    
    /* Width and height */
    int w_;
    int h_;
    
    /* Grid cell size and fluid density */
    double hx_;
    double density_;
    
    /* Arrays for: */
    VectorXd r_; /* Residual Vector */
    VectorXd p_; /* Pressure solution */
    VectorXd z_; /* Auxiliary Vector */
    VectorXd s_; /* Search Vector */
    VectorXd precon_; /* Preconditioner */

    /* Pressure Matrix stored as 3 vectors */
    VectorXd Adiag_;
    VectorXd Aplusi_;
    VectorXd Aplusj_;
    
    /* List of solid bodies to consider in the simulation */
    bool dirty;
    vector<SolidBody *> bodies_;
    
    /* Builds the pressure right hand side as the negative divergence */
    void buildResidual() {
        auto &cell = d_->cell();

        for (int y = 0, idx = 0; y < h_; y++) {
            for (int x = 0; x < w_; x++, idx++) {
                if (cell(idx) == CELL_FLUID) {
                    r_(idx) = (u_->at(x + 1, y) - u_->at(x, y) +
                               v_->at(x, y + 1) - v_->at(x, y)) / hx_ * -1;
                } else
                    r_(idx) = 0;
            }
        }
    }
    
    /* Builds the pressure matrix. Since the matrix is very sparse and
     * symmetric, it allows for memory friendly storage.
     */
    void buildA(double dt) {
        double scale = dt / (density_ * hx_ * hx_);
        auto &cell = d_->cell();
        Adiag_.setZero();
        Aplusi_.setZero();
        Aplusj_.setZero();

        for (int y = 0; y < h_; y++) {
            for (int x = 0; x < w_; x++) {
                int idx = x + y * w_;

                if (cell(idx) != CELL_FLUID)
                    continue;

                if (x < w_ - 1 && cell(idx + 1) == CELL_FLUID) {
                    Adiag_ (idx    ) +=  scale;
                    Adiag_ (idx + 1) +=  scale;
                    Aplusi_(idx    )  = -scale;
                }
                if (y < h_ - 1 && cell(idx + w_) == CELL_FLUID) {
                    Adiag_ (idx     ) +=  scale;
                    Adiag_ (idx + w_) +=  scale;
                    Aplusj_(idx     )  = -scale;
                }
            }
        }
    }
    
    /* Builds the modified incomplete Cholesky preconditioner */
    void buildPrecon() {
        const double tau = 0.97;
        const double sigma = 0.25;
        auto &cell = d_->cell();

        for (int y = 0; y < h_; y++) {
            for (int x = 0; x < w_; x++) {
                int idx = x + y * w_;

                if (cell(idx) != CELL_FLUID)
                    continue;

                double e = Adiag_(idx);

                if (x > 0 && cell(idx - 1) == CELL_FLUID) {
                    double px = Aplusi_(idx - 1) * precon_(idx - 1);
                    double py = Aplusj_(idx - 1) * precon_(idx - 1);
                    e = e - (px * px + tau * px * py);
                }

                if (y > 0 && cell(idx - w_) == CELL_FLUID) {
                    double px = Aplusi_(idx - w_) * precon_(idx - w_);
                    double py = Aplusj_(idx - w_) * precon_(idx - w_);
                    e = e - (py * py + tau * px * py);
                }

                // Not sure why add this guard?
                if (e < sigma * Adiag_(idx))
                    e = Adiag_(idx);

                precon_(idx) = 1.0 / sqrt(e + 1.0e-30);
            }
        }
    }
    
    /* Apply preconditioner to vector `a' and store it in `dst' */
    void applyPrecon(VectorXd &dst, VectorXd &a) {
        auto &cell = d_->cell();

        for (int y = 0; y < h_; y++) {
            for (int x = 0; x < w_; x++) {
                int idx = x + y * w_;

                if (cell(idx) != CELL_FLUID)
                    continue;

                double t = a(idx);

                if (x > 0 && cell(idx - 1) == CELL_FLUID)
                    t -= Aplusi_(idx -  1) * precon_(idx -  1) * dst(idx -  1);
                if (y > 0 && cell(idx - w_) == CELL_FLUID)
                    t -= Aplusj_(idx - w_) * precon_(idx - w_) * dst(idx - w_);

                dst(idx) = t * precon_(idx);
            }
        }

        for (int y = h_ - 1; y >= 0; y--) {
            for (int x = w_ - 1; x >= 0; x--) {
                int idx = x + y * w_;

                if (cell(idx) != CELL_FLUID)
                    continue;

                double t = dst(idx);

                if (x < w_ - 1 && cell(idx + 1) == CELL_FLUID)
                    t -= Aplusi_(idx) * precon_(idx) * dst(idx +  1);
                if (y < h_ - 1 && cell(idx + w_) == CELL_FLUID)
                    t -= Aplusj_(idx) * precon_(idx) * dst(idx + w_);

                dst(idx) = t * precon_(idx);
            }
        }
    }
    
    /* Multiplies internal pressure matrix with vector `b' and stores the result in `dst' */
    void applyA(VectorXd &dst, VectorXd &b) {
        for (int y = 0; y < h_; y++) {
            for (int x = 0; x < w_; x++) {
                int idx = x + y * w_;

                double t = Adiag_(idx) * b(idx);
                
                if (x > 0)
                    t += Aplusi_(idx - 1) * b(idx - 1);
                if (y > 0)
                    t += Aplusj_(idx - w_) * b(idx - w_);
                if (x < w_ - 1)
                    t += Aplusi_(idx) * b(idx + 1);
                if (y < h_ - 1)
                    t += Aplusj_(idx) * b(idx + w_);

                dst(idx) = t;
            }
        }
    }
    
    /* Performs the pressure solve using PCG.
     * The solver will run as long as it takes to get the relative error below
     * a threshold, but will never exceed `limit' iterations
     */
    void project(int limit, double dt, double tolerance) {
        buildA(dt);
        buildPrecon();
        // std::cout << "precon_ = [\n" << precon_ <<"]\n";
        // std::cout << "Adiag_ = [\n" << Adiag_ <<"]\n";
        // std::cout << "Aplusi_ = [\n" << Aplusi_ <<"]\n";
        // std::cout << "Aplusj_ = [\n" << Aplusj_ <<"]\n";
        
        p_.setZero();
        buildResidual();
        // std::cout << "r_ = [\n" << r_ <<"]\n";

        if (r_.norm() < tolerance) {
            return;
        }

        // z_.setZero();
        applyPrecon(z_, r_);

        s_ = z_;
        double sigma = z_.dot(r_);

        if (sigma == 0 || isnan(sigma)) {
            return;
        }
        // std::cout << "sigma = " << sigma << "\n";

        for (int iter = 0; iter < limit; iter++) {
            applyA(z_, s_);
            // std::cout << "z_ = [\n" << z_ <<"]\n";
            // std::cout << "s_ = [\n" << s_ <<"]\n";
            double zs = z_.dot(s_);
            // std::cout << "zs = " << zs << "\n";
            if (zs == 0 || isnan(zs)) {
                break;
            }
                
            double alpha = sigma / zs;

            // std::cout << "alpha = " << alpha << "\n";
            if (alpha == 0) {
                break;
            }

            p_ += alpha * s_;
            // std::cout << "p_ = [\n" << p_ <<"\n";
            r_ -= alpha * z_;

            if (r_.norm() < tolerance) {
                break;
            }

            applyPrecon(z_, r_);
            double sigma_ = z_.dot(r_);
            s_ = z_ + sigma_ / sigma * s_;
            sigma = sigma_;
        }

    }
    
    /* Applies the computed pressure to the velocity field */
    void applyPressure(double dt) {
        double scale = dt / (density_ * hx_);
        auto &cell = d_->cell();
        
        for (int y = 0; y < h_; y++) {
            for (int x = 0; x < w_; x++) {
                int idx = x + y * w_;

                if (cell(idx) != CELL_FLUID)
                    continue;
                u_->at(x,     y    ) -= scale * p_(x + y * w_);
                u_->at(x + 1, y    ) += scale * p_(x + y * w_);
                v_->at(x,     y    ) -= scale * p_(x + y * w_);
                v_->at(x,     y + 1) += scale * p_(x + y * w_);
            }
        }
        
        for (int y = 0; y < h_; y++)
            u_->at(0, y) = u_->at(w_, y) = 0.0;

        for (int x = 0; x < w_; x++)
            v_->at(x, 0) = v_->at(x, h_) = 0.0;
    }

    
    /* Sets all velocity cells bordering solid cells to the solid velocity */
    void setBoundaryCondition() {
        auto &cell = d_->cell();
        auto &body = d_->body();
        
        for (int y = 0, idx = 0; y < h_; y++) {
            for (int x = 0; x < w_; x++, idx++) {
                if (cell(idx) == CELL_SOLID) {
                    const SolidBody &b = *bodies_[body(idx)];
                    
                    u_->at(x, y) = b.velocityX(x*hx_, (y + 0.5)*hx_);
                    v_->at(x, y) = b.velocityY((x + 0.5)*hx_, y*hx_);
                    u_->at(x + 1, y) = b.velocityX((x + 1.0)*hx_, (y + 0.5)*hx_);
                    v_->at(x, y + 1) = b.velocityY((x + 0.5)*hx_, (y + 1.0)*hx_);
                }
            }
        }
        
        for (int y = 0; y < h_; y++)
            u_->at(0, y) = u_->at(w_, y) = 0.0;

        for (int x = 0; x < w_; x++)
            v_->at(x, 0) = v_->at(x, h_) = 0.0;
    }
    
public:
    FluidSolver(int w, int h, double density, double hx) : w_(w), h_(h), density_(density), hx_(hx) {
        d_ = new MACGrid(w_,     h_,     0.5, 0.5, hx_);
        u_ = new MACGrid(w_ + 1, h_,     0.0, 0.5, hx_);
        v_ = new MACGrid(w_,     h_ + 1, 0.5, 0.0, hx_);
        
        r_.resize(w_ * h_);
        r_.setZero();
        p_.resize(w_ * h_);
        p_.setZero();
        z_.resize(w_ * h_);
        z_.setZero();
        s_.resize(w_ * h_);
        s_.setZero();
        precon_.resize(w_ * h_);
        precon_.setZero();

        Adiag_.resize(w_ * h_);
        Adiag_.setZero();
        Aplusi_.resize(w_ * h_);
        Aplusi_.setZero();
        Aplusj_.resize(w_ * h_);
        Aplusj_.setZero();

        dirty = true;
    }
    
    ~FluidSolver() {
        delete d_;
        delete u_;
        delete v_;
    }
    
    void update(double dt, int maxIter, double maxError, BodyForces &f) {

        if (dirty) {

            d_->fillSolidFields(bodies_);
            u_->fillSolidFields(bodies_);
            v_->fillSolidFields(bodies_);
            
            setBoundaryCondition();

            dirty = false;

        }

        
        d_->extrapolate();
        u_->extrapolate();
        v_->extrapolate();

        setBoundaryCondition();


        d_->advect(dt, *u_, *v_, bodies_);
        u_->advect(dt, *u_, *v_, bodies_);
        v_->advect(dt, *u_, *v_, bodies_);
        
        /* Make effect of advection visible, since it's not an in-place operation */
        d_->flip();
        u_->flip();
        v_->flip();
        // std::cout << "1 iter done...\n";

        for (int y = 1; y < h_ - 1; y++) {
            for (int x = 0; x < w_; x++) {
                int idx = x + y * w_;
                v_->at(x, y) += f.g * dt;
            }
        }

        project(maxIter, dt, maxError);
        // std::cout << "p_ = " << p_ <<"\n";
        applyPressure(dt);

        // for (auto b : bodies_) {
        //     b->update(dt);
        // }
    }
    
    /* Set density and x/y velocity in given rectangle to d/u/v, respectively */
    void addInflow(double x, double y, double w, double h, double d, double u, double v) {
        d_->addInflow(x, y, x + w, y + h, d);
        u_->addInflow(x, y, x + w, y + h, u);
        v_->addInflow(x, y, x + w, y + h, v);
    }
    
    double atImage(int x, int y) {
        auto &cell = d_->cell();
        if (cell(x + y * w_) == CELL_SOLID)
            return -1.0;
        return d_->src()[y * w_ + x];
    }

    void addBody(SolidBody *b) {
        bodies_.push_back(b);
        dirty = true;
    }
};

