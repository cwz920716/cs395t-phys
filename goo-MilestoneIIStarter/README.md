# Name
Wenzhi Cui

# World of Goo Milestone II Start Code

A simple simulation framework using libigl and cmake. Based on Alec Jacobson's libigl example project. This project contains some boilerplate that sets up a physical simulation to run in its own thread, with rendering provided by libigl, as well as a partial implementation of Milestone I with the features relevant to Milestone II.

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

A glfw app should launch displaying the goo world and GUI.

## Dependencies

The only dependencies are stl, eigen, [libigl](libigl.github.io/libigl/) and
the dependencies of the `igl::viewer::Viewer` (mandatory: glfw and
opengl, optional: nanogui and nanovg).

We recommend you to install libigl using git via:

    git clone --recursive https://github.com/libigl/libigl.git

## Supported Feature.
I have only implemented all the required features.

## Known Issues
1. The flexible rods simulation will explode when there are two nodes are very close to each other or there is a flexible triangle which is almost "flat", i.e., an angle close to 180 degree. The explosion happens no matter what numerical method we use.

2. When two nodes connected by flexible rods are very close and one of them is fixed, then the fixed one might disappear during some experiments.

3. When the number of segments and Bending stiffness is too large, flexible rods will also blow up.

