# Name
Wenzhi Cui

# World of Goo Milestone I Start Code

A simple simulation framework using libigl and cmake. Based on Alec Jacobson's libigl example project. This project contains some boilerplate that sets up a physical simulation to run in its own thread, with rendering provided by libigl.

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example_bin` binary.

## Run

From within the `build` directory just issue:

    ./goo_bin

A glfw app should launch displaying a 3D cube.

## Dependencies

The only dependencies are stl, eigen, [libigl](libigl.github.io/libigl/) and
the dependencies of the `igl::viewer::Viewer` (mandatory: glfw and
opengl, optional: nanogui and nanovg).

We recommend you to install libigl using git via:

    git clone --recursive https://github.com/libigl/libigl.git

I installed this library in the same location as the goo-MilestonelStarter folder.

## Supported Feature.
I have only implemented all the required features. Most the required features are 
implemented in their own function, like springForceHeissan() will compute the Heissan 
of the spring potentials.

My floor force is implemented as follows: if the particle is below floor, then it will
be influenced by a force F(x, y) = (0, 5m * fabs(g)) + (-2v(x), -2v(y)) and v(x, y) = (q_{i}(x, y) - q_{i-1}(x, y)) / h.
The first part is like an reversed gravity (I also called fabs() on g to make sure this force is always
upward) and the second part is like a drag which will slow the particle down. 

Also, my program will output some lines like "sprinf snapped" or "particle left scene" when
they get deleted. It is mainly for debugging so feel free to comment out if you must.

I have tested on lab machines and included a goo_bin in the build folder in case it does not compile due to some reasons.

## Known Issues
1. I have a segment fault when I try to close the screen.

2. The simulation slows down when the number of objects increases (more than 20).
Espeicially when using the Newton method.
