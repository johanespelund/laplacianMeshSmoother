# laplacianMeshSmoother
An OpenFOAM utility to use laplacian smoothing. **WIP!**

## Installation
Tested using OpenFOAM v2206.

## Usage

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.1.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      laplacianSmoothDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Prevent moving points in the boundary normal direction, for points a distance from the boundaries.
// Not working properly. Instead, refineWallLayer is used after smoothing.
preserveBoundaryLayer #calc "$t_BL / 2"; //#eval {10*$wall_cell_size}; 
iters 20;
smoothFactor 0.8;
boundaryNormalFreq 0;
boundaryNormalPatches ("wall_.*" );

constrainedPoints
(
  {
    type patch;
    patch wall;
    constraintType fixed;
  }

  {
    type patch;
    patch wedgePos;
    constraintType constDir;
    direction (-0.00872654  0  0.999962);
  }
  {
    type patch;
    patch wedgeNeg;
    constraintType constDir;
    direction (-0.00872654  0  -0.999962);
  }
);
```

## OpenFOAM version notes
- OpenFOAM.com and OpenFOAM.org use different APIs in some places.
- This codebase includes conditional handling for both flavors.
- `laplacianMeshSmoother.C` currently has `#define ORG_VERSION` enabled by default.

Clone repo and run `wmake`. The executable will be located in `$FOAM_USER_APPBIN` and can be called with `laplacianMeshSmoother` when OpenFOAM is loaded.

TODO:
- [ ] Add automatic selection of points based on non-orthogonality and/or skewness
- [ ] Fix distortion of boundary layer
  - Add constraint that limits movement in the boundary normal direction, based on the distance from the boundary. This could eliminate the need for the yconst and rconst constraints.
  - Another option is to just constratin the boundary points and let other points be free. This would at least eliminate the need for hardcoded geometry constraints such as yconst and rconst.
    If points that belong to different patches are marked as fixedPoints, the boundary should not change.
  - Current strategy is to use hardcoded constraints and fixedPoints, and then use refineWallLayer if the first cell height has increase significantly.
