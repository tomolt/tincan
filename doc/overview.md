An Overview of the internals of tincan
======================================

* TBD Broadphase
  - persistent, incremental Sweep-and-prune
  - Based on insertion-sort
  - Largely based on report by Terdiman
* Narrowphase / Collision detection
  - Minkowski Portal Refinement
  - TBD SIMD Polytope support function
  - TBD Warmstarting with old portal
  - Contact point is also found via MPR
  - TBD Multithreading
* Contact manifold construction
  - Incrementally build Manifold via contact caching
  - Contact point of portal is found via barycentric coordinates
  - TBD SIMD contact point calculation
* TBD Simulation island detection
  - Depth-first-search over collision graph
  - Sleeping?
* Constraint resolution
  - Sequential Impulses
  - Gauss-Seidel solver
  - Precompute effective mass
  - Rematerialize Jacobian when needed
  - TBD Warmstarting of Lagrange magnitudes (needs projection onto new Jacobian)
  - TBD Process islands on separate threads
  - TBD SIMD? Perhaps Gauss-Seidel-Jacobi hybrid
  - TBD transform inverse Inertia etc. in world-space beforehand
* Position integration
  - Semi-implicit Euler integration

