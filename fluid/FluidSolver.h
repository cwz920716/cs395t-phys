#pragma once

#include <algorithm>
#include <math.h>
#include <stdio.h>

using namespace std;

class FluidQuantity {
    /* Memory buffers for fluid quantity */
    double *_src;
    double *_dst;

    /* Width and height */
    int _w;
    int _h;
    /* X and Y offset from top left grid cell.
     * This is (0.5,0.5) for centered quantities such as density,
     * and (0.0, 0.5) or (0.5, 0.0) for jittered quantities like the velocity.
     */
    double _ox;
    double _oy;
    /* Grid cell size */
    double _hx;
    
    /* Linear intERPolate between a and b for x ranging from 0 to 1 */
    double lerp(double a, double b, double x) const {
        return a*(1.0 - x) + b*x;
    }
    
    /* Simple forward Euler method for velocity integration in time */
    void euler(double &x, double &y, double timestep, const FluidQuantity &u, const FluidQuantity &v) const {
        double uVel = u.lerp(x, y)/_hx;
        double vVel = v.lerp(x, y)/_hx;
        
        x -= uVel*timestep;
        y -= vVel*timestep;
    }
    
public:
    FluidQuantity(int w, int h, double ox, double oy, double hx)
            : _w(w), _h(h), _ox(ox), _oy(oy), _hx(hx) {
        _src = new double[_w*_h];
        _dst = new double[_w*_h];
                
        memset(_src, 0, _w*_h*sizeof(double));
    }
    
    ~FluidQuantity() {
        delete[] _src;
        delete[] _dst;
    }
    
    void flip() {
        swap(_src, _dst);
    }
    
    const double *src() const {
        return _src;
    }
    
    /* Read-only and read-write access to grid cells */
    double at(int x, int y) const {
        return _src[x + y*_w];
    }
    
    double &at(int x, int y) {
        return _src[x + y*_w];
    }
    
    /* Linear intERPolate on grid at coordinates (x, y).
     * Coordinates will be clamped to lie in simulation domain
     */
    double lerp(double x, double y) const {
        x = min(max(x - _ox, 0.0), _w - 1.001);
        y = min(max(y - _oy, 0.0), _h - 1.001);
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
        for (int iy = 0, idx = 0; iy < _h; iy++) {
            for (int ix = 0; ix < _w; ix++, idx++) {
                double x = ix + _ox;
                double y = iy + _oy;
                
                /* First component: Integrate in time */
                euler(x, y, timestep, u, v);
                
                /* Second component: Interpolate from grid */
                _dst[idx] = lerp(x, y);
            }
        }
    }
    
    /* Sets fluid quantity inside the given rect to value `v' */
    void addInflow(double x0, double y0, double x1, double y1, double v) {
        int ix0 = (int)(x0/_hx - _ox);
        int iy0 = (int)(y0/_hx - _oy);
        int ix1 = (int)(x1/_hx - _ox);
        int iy1 = (int)(y1/_hx - _oy);
        
        for (int y = max(iy0, 0); y < min(iy1, _h); y++)
            for (int x = max(ix0, 0); x < min(ix1, _h); x++)
                if (fabs(_src[x + y*_w]) < fabs(v))
                    _src[x + y*_w] = v;
    }
};

/* Fluid solver class. Sets up the fluid quantities, forces incompressibility
 * performs advection and adds inflows.
 */
class FluidSolver {
    /* Fluid quantities */
    FluidQuantity *_d;
    FluidQuantity *_u;
    FluidQuantity *_v;
    
    /* Width and height */
    int _w;
    int _h;
    
    /* Grid cell size and fluid density */
    double _hx;
    double _density;
    
    /* Arrays for: */
    double *_r; /* Right hand side of pressure solve */
    double *_p; /* Pressure solution */
    
    
    /* Builds the pressure right hand side as the negative divergence */
    void buildRhs() {
        double scale = 1.0/_hx;
        
        for (int y = 0, idx = 0; y < _h; y++) {
            for (int x = 0; x < _w; x++, idx++) {
                _r[idx] = -scale*(_u->at(x + 1, y) - _u->at(x, y) +
                                  _v->at(x, y + 1) - _v->at(x, y));
            }
        }
    }
    
    /* Performs the pressure solve using Gauss-Seidel.
     * The solver will run as long as it takes to get the relative error below
     * a threshold, but will never exceed `limit' iterations
     */
    void project(int limit, double timestep) {
        double scale = timestep/(_density*_hx*_hx);
        
        double maxDelta;
        for (int iter = 0; iter < limit; iter++) {
            maxDelta = 0.0;
            for (int y = 0, idx = 0; y < _h; y++) {
                for (int x = 0; x < _w; x++, idx++) {
                    int idx = x + y*_w;
                    
                    double diag = 0.0, offDiag = 0.0;
                    
                    /* Here we build the matrix implicitly as the five-point
                     * stencil. Grid borders are assumed to be solid, i.e.
                     * there is no fluid outside the simulation domain.
                     */
                    if (x > 0) {
                        diag    += scale;
                        offDiag -= scale*_p[idx - 1];
                    }
                    if (y > 0) {
                        diag    += scale;
                        offDiag -= scale*_p[idx - _w];
                    }
                    if (x < _w - 1) {
                        diag    += scale;
                        offDiag -= scale*_p[idx + 1];
                    }
                    if (y < _h - 1) {
                        diag    += scale;
                        offDiag -= scale*_p[idx + _w];
                    }

                    double newP = (_r[idx] - offDiag)/diag;
                    
                    maxDelta = max(maxDelta, fabs(_p[idx] - newP));
                    
                    _p[idx] = newP;
                }
            }

            if (maxDelta < 1e-5) {
                return;
            }
        }
    }
    
    /* Applies the computed pressure to the velocity field */
    void applyPressure(double timestep) {
        double scale = timestep/(_density*_hx);
        
        for (int y = 0, idx = 0; y < _h; y++) {
            for (int x = 0; x < _w; x++, idx++) {
                _u->at(x,     y    ) -= scale*_p[idx];
                _u->at(x + 1, y    ) += scale*_p[idx];
                _v->at(x,     y    ) -= scale*_p[idx];
                _v->at(x,     y + 1) += scale*_p[idx];
            }
        }
        
        for (int y = 0; y < _h; y++)
            _u->at(0, y) = _u->at(_w, y) = 0.0;
        for (int x = 0; x < _w; x++)
            _v->at(x, 0) = _v->at(x, _h) = 0.0;
    }
    
public:
    FluidSolver(int w, int h, double density, double hx) : _w(w), _h(h), _density(density), _hx(hx) {
        _d = new FluidQuantity(_w,     _h,     0.5, 0.5, _hx);
        _u = new FluidQuantity(_w + 1, _h,     0.0, 0.5, _hx);
        _v = new FluidQuantity(_w,     _h + 1, 0.5, 0.0, _hx);
        
        _r = new double[_w*_h];
        _p = new double[_w*_h];
        
        memset(_p, 0, _w*_h*sizeof(double));
    }
    
    ~FluidSolver() {
        delete _d;
        delete _u;
        delete _v;
        
        delete[] _r;
        delete[] _p;
    }
    
    void update(double timestep) {
        buildRhs();
        project(600, timestep);
        applyPressure(timestep);
        
        _d->advect(timestep, *_u, *_v);
        _u->advect(timestep, *_u, *_v);
        _v->advect(timestep, *_u, *_v);
        
        /* Make effect of advection visible, since it's not an in-place operation */
        _d->flip();
        _u->flip();
        _v->flip();
    }
    
    /* Set density and x/y velocity in given rectangle to d/u/v, respectively */
    void addInflow(double x, double y, double w, double h, double d, double u, double v) {
        _d->addInflow(x, y, x + w, y + h, d);
        _u->addInflow(x, y, x + w, y + h, u);
        _v->addInflow(x, y, x + w, y + h, v);
    }
    
    /* Returns the maximum allowed timestep. Note that the actual timestep
     * taken should usually be much below this to ensure accurate
     * simulation - just never above.
     */
    double maxTimestep() {
        double maxVelocity = 0.0;
        for (int y = 0; y < _h; y++) {
            for (int x = 0; x < _w; x++) {
                /* Average velocity at grid cell center */
                double u = _u->lerp(x + 0.5, y + 0.5);
                double v = _v->lerp(x + 0.5, y + 0.5);
                
                double velocity = sqrt(u*u + v*v);
                maxVelocity = max(maxVelocity, velocity);
            }
        }
        
        /* Fluid should not flow more than two grid cells per iteration */
        double maxTimestep = 2.0*_hx/maxVelocity;
        
        /* Clamp to sensible maximum value in case of very small velocities */
        return min(maxTimestep, 1.0);
    }
    
    /* Convert fluid density to RGBA image */
    double atImage(int x, int y) {
        return _d->src()[y * _w + x];
    }
};

