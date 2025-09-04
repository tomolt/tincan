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

Intersection tests between convex polytopes are performed through the XenoCollide algorithm, aka Minkowski Portal Refinement.
The contact normal is found through iteration of the MPR algorithm.
Antiparallel faces along the contact normal are chosen as contact features.
Every frame, a new contact manifold is built through projection and Sutherland-Hogdman clipping of the contact features.
Constraints are resolved by a sequential impulses solver.

