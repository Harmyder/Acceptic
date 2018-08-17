The application based on [On the Separation of a Polyhedron from Its Single-Part Mold][1].

It implements 2D LP-solver based on Seidel’s randomized incremental algorithm which runs in expected linear time.

To find a cover for S2 by faces' hemispheres, it needs to find a cover for each of four predefined hemispheres which cover S2. 

For every such predefined hemisphere the task boils down to find a half planes coverage in the projection plane `z = 1`. The algorithms proceeds with the following two steps:

1. 
   1. Project all upward-looking faces onto the projection plane.
   2. Choose arbitrary objective vector `c`. Run LP two times. First with `c`, second with `-c`. The solver returns not only the best value but tight halfplanes as well. There are could 0, 1 or 2 of them. If after these two runs i have coverage from these half planes then step 2 is not needed, otherwise continue with step 2.
   
2. By now it has either one or two halfplanes which it can augment to have a coverage. It needs to augment with half planes formed by downward-looking faces. For this end it projects them on the same projection plane. Next it have two options:
   1. Find one-augment. Easy, just find a half plane which has normal negative to one of the half planes normals it has so far. And check if it forms a coverage. If so, it's done, otherwise proceed with 2.2.
   2. Find two-augment. That's more complex, but still algorithm is similar to what is done in the LP solver. You can look at ASCII-graphics in `Casting.cpp` for `FindTwoMoreForCover`.
   
After the algorithm has all candidates in place it uses function `CheckCandidates` to find a good pullout for every candidate. Here it again resorts to LP with `c` and `-c`. Every run of the solver returns border point which allows to pull out the form. The algorithm returns average of the two in the hope that it is good enough.


PS. If `Solver::Randomize` is called then solver runs in expected linear time, and so the whole algorithm. My goal was to keep algorithm theoretical complexity in these bounds.


[1]: https://arxiv.org/pdf/1708.04203.pdf
