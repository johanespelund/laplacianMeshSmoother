# laplacianMeshSmoother
An OpenFOAM utility to use laplacian smoothing. WIP!

TODO:
- [ ] Add automatic selection of points based on non-orthogonality and/or skewness
- [ ] Fix distortion of boundary layer
  - Add constraint that limits movement in the boundary normal direction, based on the distance from the boundary. This could eliminate the need for the yconst and rconst constraints.
  - Another option is to just constratin the boundary points and let other points be free. This would at least eliminate the need for hardcoded geometry constraints such as yconst and rconst.
    If points that belong to different patches are marked as fixedPoints, the boundary should not change.
