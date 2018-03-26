# Name

Wenzhi Cui (WC8348)

# SevenEves Milestone I Start Code

A simple simulation framework using libigl and cmake. Based on Alec Jacobson's libigl example project. This project contains some boilerplate that sets up a physical simulation to run in its own thread, with rendering provided by libigl.

## Compile
I have modified CMakeLists.txt but it should work the same way.

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `example_bin` binary.

## Run

From within the `build` directory just issue:

    ./seven_bin

A glfw app should launch displaying a GUI. The code will try to load scene files from either
the ./scenes or ../scenes folder, so you need to run the binary from either the project root
folder, or a build subdirectory.

## Supported Feature.
I have only implemented all the required features.

## Dependencies

The only dependencies are stl, eigen, [libigl](libigl.github.io/libigl/) and
the dependencies of the `igl::viewer::Viewer` (glfw and opengl).

We recommend you to install libigl using git via:

    git clone --recursive https://github.com/libigl/libigl.git
