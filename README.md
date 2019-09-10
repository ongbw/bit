bit
===

C++ Boundary Integral Treecode

Boundary Integral Treecode (BIT) methods apply fast summation methods to aproximate the boundary integral formulation of potential solutions.  Specifically, a treecode alogrithm is employed to reduce the operation count of evaluating the summation (arising from the integral solution) from O(N) to O(log N).  Below, we describe the algorithm applied to canonical example of charged particles.  The algorithm can be more widely applied to different problems.

The treecode algorithm divides the particles into a hierarchy of clusters, and the particle-particle interactions are placed with particle-cluster interactions which are approximated using multi-pole expansions.  Barnes and Hut (1986) used mono-pole approximations and a divide-and-conquer evaluation strategy, while Greengard and Rokhlin (1987) used higher order spherical harmonics expansions and a more sophisticated evaluation procedure.  Treecode algorithms have been very successful in particle simulations and there is ongoing interest in optimizing their performance.

Early version of this treecode software implemented the treecode as described in
<http://www.math.lsa.umich.edu/~krasny/paper_jcp_2001.pdf>

The latest version of this treecode implements a two-parameter family of regularized kernels described in
<http://dx.doi.org/10.1007/s10915-016-0336-0>
