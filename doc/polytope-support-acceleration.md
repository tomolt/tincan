Polytope Support Function Acceleration Datastructure
===

I require an efficient way of querying the farthest point within a polytope along a given direction vector.
Additionally, I want to be able to query all points in the polytope that lie on a given plane on the hull of the polytope.

The core idea is that I can precompute the farthest points along a number of directions.
Then, when answering a query, I can look up a precomputed farthest point for a similar direction,
project that point onto the search direction, and use that distance as a lower bound for the search.

One can take this idea a step further by precomputing points for a 'grid' of directions, e.g. on the surface of the unit sphere.
Then any query direction falls within a square of the unit sphere, where each corner of the square has a precomputed point.
I call these squares 'buckets'.
For each bucket, I want to precompute the list of all points of the polytope which may be further on any direction covered by the square than its corners.
From there, answering queries will be easy.
One notable advantage of this scheme is that polytope points which do not lie on the hull wont bloat the datastructure.
Instead of a plain grid, I might want to upgrade to a quad-tree.
That would allow me to subdivide buckets until all lists conform to a fixed maximum length, e.g. 4 points.
The natural fanout of 4 together with the maximum length lend themselves well to perform some easy 4-way SIMD computations.

In practice, I will use a different direction grid that is easier to compute.
I will project direction vectors onto the surface of an octohedron, and then unfold them into a square.

But how can we precompute the bucket lists?
First, an observation in the 2D case, where buckets are lines between two corners:
If both corners happen to have the same point, then 
