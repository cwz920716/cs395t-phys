#pragma once

#include <algorithm>
#include <math.h>
#include <stdio.h>

using namespace std;
using namespace Eigen;

class FluidQuantity {
    /* Memory buffers for fluid quantity */
    double *src_;
    double *dst_;

    /* Width and height */
    int w_;
    int h_;
    /* X and Y offset from top left grid cell.
     * This is (0.5,0.5) for centered quantities such as density,
     * and (0.0, 0.5) or (0.5, 0.0) for jittered quantities like the velocity.
     */
    double ox_;
    double _oy;
    /* Grid cell size */
    double hx_;
    
    /* Linear intERPolate between a and b for x ranging from 0 to 1 */
    double lerp(double a, double b, double x) const {
        return a*(1.0 - x) + b*x;
    }
    
    /* Simple forward Euler method for velocity integration in time */
    void euler(double &x, double &y, double timestep, const FluidQuantity &u, const FluidQuantity &v) const {
        double uVel = u.lerp(x, y)/hx_;
        double vVel = v.lerp(x, y)/hx_;
        
        x -= uVel*timestep;
        y -= vVel*timestep;
    }
    
public:
    FluidQuantity(int w, int h, double ox, double oy, double hx)
            : w_(w), h_(h), ox_(ox), _oy(oy), hx_(hx) {
        src_ = new double[w_*h_];
        dst_ = new double[w_*h_];
                
        memset(src_, 0, w_*h_*sizeof(double));
    }
    
    ~FluidQuantity() {
        delete[] src_;
        delete[] dst_;
    }
    
    void flip() {
        swap(src_, dst_);
    }
    
    const double *src() const {
        return src_;
    }
    
    /* Read-only and read-write access to grid cells */
    double at(int x, int y) const {
        return src_[x + y*w_];
    }
    
    double &at(int x, int y) {
        return src_[x + y*w_];
    }
    
    /* Linear intERPolate on grid at coordinates (x, y).
     * Coordinates will be clamped to lie in simulation domain
     */
    double lerp(double x, double y) const {
        x = min(max(x - ox_, 0.0), w_ - 1.001);
        y = min(max(y - _oy, 0.0), h_ - 1.001);
        int ix = (int)x;
        int iy = (int)y;
        x -= ix;
        y -= iy;
        
        double x00 = at(ix + 0, iy + 0), x10 = at(ix + 1, iy + 0);
        double x01 = at(ix + 0, iy + 1), x11 = at(ix + 1, iy + 1);
        
        return lerp(lerp(x00, x10, x), lerp(x01, x11, x), y);
    }
    
    /* Advect grid in velocity field u, v with given timestep */
    void advect(double timestep, const FluidQuantity &u, const FluidQuantity &v) {
        for (int iy = 0, idx = 0; iy < h_; iy++) {
            for (int ix = 0; ix < w_; ix++, idx++) {
                double x = ix + ox_;
                double y = iy + _oy;
                
                /* First component: Integrate in time */
                euler(x, y, timestep, u, v);
                
                /* Second component: Interpolate from grid */
                dst_[idx] = lerp(x, y);
            }
        }
    }
    
    /* Sets fluid quantity inside the given rect to value `v' */
    void addInflow(double x0, double y0, double x1, double y1, double v) {
        int ix0 = (int)(x0/hx_ - ox_);
        int iy0 = (int)(y0/hx_ - _oy);
        int ix1 = (int)(x1/hx_ - ox_);
        int iy1 = (int)(y1/hx_ - _oy);
        
        for (int y = max(iy0, 0); y < min(iy1, h_); y++)
            for (int x = max(ix0, 0); x < min(ix1, h_); x++)
                if (fabs(src_[x + y*w_]) < fabs(v))
                    src_[x + y*w_] = v;
    }
};

/* Fluid solver class. Sets up the fluid quantities, forces incompressibility
 * performs advection and adds inflows.
 */
class FluidSolver {
    /* Fluid quantities */
    FluidQuantity *d_;
    FluidQuantity *u_;
    FluidQuantity *v_;
    
    /* Width and height */
    int w_;
    int h_;
    
    /* Grid cell size and fluid density */
    double hx_;
    double density_;
    
    /* Arrays for: */
    double *r_; /* Right hand side of pressure solve */
    double *p_; /* Pressure solution */
    
    
    /* Builds the pressure right hand side as the negative divergence */
    void buildRhs() {
        double scale = 1.0/hx_;
        
        for (int y = 0, idx = 0; y < h_; y++) {
            for (int x = 0; x < w_; x++, idx++) {
                r_[idx] = -scale*(u_->at(x + 1, y) - u_->at(x, y) +
                                  v_->at(x, y + 1) - v_->at(x, y));
            }
        }
    }
    
    /* Performs the pressure solve using Gauss-Seidel.
     * The solver will run as long as it takes to get the relative error below
     * a threshold, but will never exceed `limit' iterations
     */
    void project(int limit, double timestep) {
        double scale = timestep/(density_*hx_*hx_);
        
        double maxDelta;
        for (int iter = 0; iter < limit; iter++) {
            maxDelta = 0.0;
            for (int y = 0, idx = 0; y < h_; y++) {
                for (int x = 0; x < w_; x++, idx++) {
                    int idx = x + y*w_;
                    
                    double diag = 0.0, offDiag = 0.0;
                    
                    /* Here we build the matrix implicitly as the five-point
                     * stencil. Grid borders are assumed to be solid, i.e.
                     * there is no fluid outside the simulation domain.
                     */
                    if (x > 0) {
                        diag    += scale;
                        offDiag -= scale*p_[idx - 1];
                    }
                    if (y > 0) {
                        diag    += scale;
                        offDiag -= scale*p_[idx - w_];
                    }
                    if (x < w_ - 1) {
                        diag    += scale;
                        offDiag -= scale*p_[idx + 1];
                    }
                    if (y < h_ - 1) {
                        diag    += scale;
                        offDiag -= scale*p_[idx + w_];
                    }

                    double newP = (r_[idx] - offDiag)/diag;
                    
                    maxDelta = max(maxDelta, fabs(p_[idx] - newP));
                    
                    p_[idx] = newP;
                }
            }

            if (maxDelta < 1e-5) {
                return;
            }
        }
    }
    
    /* Applies the computed pressure to the velocity field */
    void applyPressure(double timestep) {
        double scale = timestep/(density_*hx_);
        
        for (int y = 0, idx = 0; y < h_; y++) {
            for (int x = 0; x < w_; x++, idx++) {
                u_->at(x,     y    ) -= scale*p_[idx];
                u_->at(x + 1, y    ) += scale*p_[idx];
                v_->at(x,     y    ) -= scale*p_[idx];
                v_->at(x,     y + 1) += scale*p_[idx];
            }
        }
        
        for (int y = 0; y < h_; y++)
            u_->at(0, y) = u_->at(w_, y) = 0.0;
        for (int x = 0; x < w_; x++)
            v_->at(x, 0) = v_->at(x, h_) = 0.0;
    }
    
public:
    FluidSolver(int w, int h, double density, double hx) : w_(w), h_(h), density_(density), hx_(hx) {
        d_ = new FluidQuantity(w_,     h_,     0.5, 0.5, hx_);
        u_ = new FluidQuantity(w_ + 1, h_,     0.0, 0.5, hx_);
        v_ = new FluidQuantity(w_,     h_ + 1, 0.5, 0.0, hx_);
        
        r_ = new double[w_*h_];
        p_ = new double[w_*h_];
        
        memset(p_, 0, w_*h_*sizeof(double));
    }
    
    ~FluidSolver() {
        delete d_;
        delete u_;
        delete v_;
        
        delete[] r_;
        delete[] p_;
    }
    
    void update(double timestep) {
        buildRhs();
        project(600, timestep);
        applyPressure(timestep);
        
        d_->advect(timestep, *u_, *v_);
        u_->advect(timestep, *u_, *v_);
        v_->advect(timestep, *u_, *v_);
        
        /* Make effect of advection visible, since it's not an in-place operation */
        d_->flip();
        u_->flip();
        v_->flip();
    }
    
    /* Set density and x/y velocity in given rectangle to d/u/v, respectively */
    void addInflow(double x, double y, double w, double h, double d, double u, double v) {
        d_->addInflow(x, y, x + w, y + h, d);
        u_->addInflow(x, y, x + w, y + h, u);
        v_->addInflow(x, y, x + w, y + h, v);
    }
    
    /* Returns the maximum allowed timestep. Note that the actual timestep
     * taken should usually be much below this to ensure accurate
     * simulation - just never above.
     */
    double maxTimestep() {
        double maxVelocity = 0.0;
        for (int y = 0; y < h_; y++) {
            for (int x = 0; x < w_; x++) {
                /* Average velocity at grid cell center */
                double u = u_->lerp(x + 0.5, y + 0.5);
                double v = v_->lerp(x + 0.5, y + 0.5);
                
                double velocity = sqrt(u*u + v*v);
                maxVelocity = max(maxVelocity, velocity);
            }
        }
        
        /* Fluid should not flow more than two grid cells per iteration */
        double maxTimestep = 2.0*hx_/maxVelocity;
        
        /* Clamp to sensible maximum value in case of very small velocities */
        return min(maxTimestep, 1.0);
    }
    
    /* Convert fluid density to RGBA image */
    double atImage(int x, int y) {
        return d_->src()[y * w_ + x];
    }
};

