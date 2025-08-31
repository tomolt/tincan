tincan physics engine
=====================

Dependencies
============

Compile-time dependencies for the core physics engine:

- A POSIX shell
- GNU make
- A GCC-compatible C11 compiler

Additionally, for the demo program:

- GLFW 3

Building
========

First, run the `configure.sh` script. This will generate a Makefile.
Then, simply run `make` to build both the physics engine and the demo program.

Features
========

- Basic 3D Rigid Body Physics
- Sliding Friction
- Convex polytope collision meshes
- Joints (for chains, ragdolls, ...)

Technology
==========

Intersection tests are performed through the XenoCollide algorithm, aka Minkowski Portal Refinement.
Collision points on the minkowski sum are found through a second MPR step.
They are converted back into world space through a roundtrip over barycentric coordinates.
Contact manifolds are cached across frames, and built up incrementally from the collision points found by MPR.

