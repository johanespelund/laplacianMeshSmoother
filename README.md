# laplacianMeshSmoother
An OpenFOAM utility to use laplacian smoothing. **WIP!**

## Installation
Tested using OpenFOAM v2206
Clone repo and run `wmake`. The executable will be located in `$FOAM_USER_APPBIN` and can be called with `laplacianMeshSmoother` when OpenFOAM is loaded.

TODO:
- [ ] Add automatic selection of points based on non-orthogonality and/or skewness
- [ ] Fix distortion of boundary layer
  - Add constraint that limits movement in the boundary normal direction, based on the distance from the boundary. This could eliminate the need for the yconst and rconst constraints.
  - Another option is to just constratin the boundary points and let other points be free. This would at least eliminate the need for hardcoded geometry constraints such as yconst and rconst.
    If points that belong to different patches are marked as fixedPoints, the boundary should not change.
  - Current strategy is to use hardcoded constraints and fixedPoints, and then use refineWallLayer if the first cell height has increase significantly.
