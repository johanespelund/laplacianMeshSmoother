/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laplacianMeshSmoother

Group
    grpMeshAdvancedUtilities

Description
    Manipulate mesh elements.

    Actions are:
        (boundary)points:
            - move

        (boundary)edges:
            - split and move introduced point

        (boundary)faces:
            - split(triangulate) and move introduced point

        edges:
            - collapse

        cells:
            - split into polygonal base pyramids around newly introduced mid
              point

    Is a bit of a loose collection of mesh change drivers.

\*---------------------------------------------------------------------------*/

#define OF12

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "Pair.H"
#include "globalIndex.H"
#include "topoSet.H"
#include "IOdictionary.H"
#include "dictionary.H"
#include "pointSet.H"
#include "faceList.H"
#include "wallDist.H"
#include "volPointInterpolation.H"
#include "wedgePolyPatch.H"
#include "wordReList.H"

#include <chrono>

#ifdef OF12

#include "systemDict.H"

#else

#include "PackedBoolList.H"


#endif

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#include "polyMesh.H"
#include "primitivePatch.H"
#include "pointField.H"

using namespace Foam;

// Function to calculate the shortest distance from a given mesh point to the nearest boundary patch
label getClosestBoundaryPoint(const polyMesh& mesh, label pointI, scalar cutoff, label patchI) {
  const vectorField& points = mesh.points();
  const point& targetPoint = points[pointI];

  scalar minDistSqr = GREAT; // Initialize with a very large number
  label closestPointIndex = -1;

  // Iterate over all boundary patches
  /* forAll(mesh.boundaryMesh(), patchI) { */
  polyPatch pp = mesh.boundaryMesh()[patchI];

  if (!pp.coupled()) { // Skip coupled patches to avoid internal faces
    labelList patchPoints = pp.meshPoints();

    // Iterate over all points in the current patch to find the closest point
    forAll(patchPoints, patchPointI) {
      scalar distSqr = magSqr(targetPoint - mesh.points()[patchPoints[patchPointI]]);

      // Update minimum distance if a closer point is found
      if (distSqr < minDistSqr) {
        minDistSqr = distSqr;
        closestPointIndex = patchPoints[patchPointI];
      }
      /* if (previousClosestPatch < patchI && distSqr < cutoffSqr) { */
      /*   /1* Info << "Previous closest patch: " << previousClosestPatch << " Current closest patch: " << patchI << endl; *1/ */
      /*   /1* Info << "Info: Multiple boundary patches within cutoff distance of point " << pointI << endl; *1/ */
      /*   return -1; */
      /* } */
      /* if (cutoff > 0.0 && distSqr < cutoffSqr) { */
      /*   belowCutoff = true; */
      /*   previousClosestPatch = patchI; */
      /* } */
    }
  }
  /* if (belowCutoff == false) */
  /* { */
  /*   previousClosestPatch++; */
  /* } */
  /* } */

  // Return the square root of the shortest squared distance found
  /* if (belowCutoff) { */
  return closestPointIndex;
  /* } */
  /* else { */
  /* return -1; */
  /* } */
}

vector getBoundaryNormal(const polyMesh& mesh, const label pointIndex) {
  /* pointField points = mesh.points(); */
  const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

  vector normal(0, 0, 0);

  // Iterate over boundary faces to find the face containing the given point index
  forAll(boundaryMesh, patchI) {
    polyPatch pp = boundaryMesh[patchI];
    forAll(boundaryMesh[patchI], faceI)
    {
      face currentFace = boundaryMesh[patchI][faceI];
      #ifdef OF12
      if (findIndex(currentFace, pointIndex) != -1)
      #else
      if (currentFace.found(pointIndex))
      #endif
      {
        vector faceNormal = boundaryMesh[patchI].faceAreas()[faceI];
        normal += faceNormal;
      }
    }
  }
  if (mag(normal) > SMALL)
  {
    return normal/mag(normal);
  }
  else {
    return vector(0, 0, 0);
  }
}


// Function to project a point onto a sphere
point projectPointOntoSphere(const vector& point, double radius) {
  // Calculate the magnitude of the point vector
  double magnitude = mag(point);

  // Normalize the point vector to get the direction
  vector direction = point / magnitude;

  // Scale the direction vector by the radius to get the coordinates of the projected point
  vector projectedPoint = direction * radius;

  return projectedPoint;
}

point projectPointOntoZCylinder(const vector& point, double radius) {
  // Project a point onto a circli in the xy-plane, keep z-coordinate
  vector projectedPoint = point;
  projectedPoint.z() = 0;
  scalar currentRadius = mag(projectedPoint);
  scalar scaleFactor = radius / currentRadius;
  projectedPoint.x() *= scaleFactor;
  projectedPoint.y() *= scaleFactor;
  projectedPoint.z() = point.z();

  return projectedPoint;
}

// Function to find neighboring points of a given point
labelList findNeighboringPoints(const polyMesh& mesh, label pointIndex) {
  labelList pointEdges = mesh.pointEdges()[pointIndex];
  labelList neighboringPointIndices;

  forAll(pointEdges, edgeI) {
    const edge& currentEdge = mesh.edges()[pointEdges[edgeI]];

    // Find the neighboring point of the current edge
    label neighboringPointIndex = (currentEdge.start() == pointIndex) ? currentEdge.end() : currentEdge.start();

    // Add the neighboring point to the list if it's not the same as the original point
    if (neighboringPointIndex != pointIndex) {
      neighboringPointIndices.append(neighboringPointIndex);
    }
  }

  // Remove duplicates
    #ifdef OF12
    Foam::sort(neighboringPointIndices);
    #else
    neighboringPointIndices = Foam::uniqueSort(neighboringPointIndices);
    #endif

  return neighboringPointIndices;
}

int main(int argc, char *argv[]) {

  argList::addNote("Boundary smoothing utility.");
  #include "addOverwriteOption.H"
  #include "addRegionOption.H"
  argList::addOption("dict", "file", "Alternative refineMeshDict");
  #include "setRootCase.H"
  #ifdef OF12
  #include "createTimeNoFunctionObjects.H"
  Foam::word meshRegionName = polyMesh::defaultRegion;
  args.optionReadIfPresent("region", meshRegionName);
  #else
  argList::noFunctionObjects(); // Never use function objects
  #include "createTime.H"
  #include "addRegionOption.H"
  #endif

  #include "createNamedPolyMesh.H"

  // Obtain dictPath here for messages
  #ifdef OF12
  fileName dictPath = args.optionLookupOrDefault<fileName>("dict", "");
  #else
  fileName dictPath = args.getOrDefault<fileName>("dict", "");
  #endif

  const word dictName("laplacianSmoothDict");
  // Dictionary to control refinement
  Info << "Reading dictionary " << dictName << " from " << dictPath << endl;

  #ifdef OF12

  const IOdictionary dict(systemDict(dictName, args, mesh));

  #else

  dictionary dict;
  #include "setSystemMeshDictionaryIO.H"
  // Read the dictionary
  dict = IOdictionary(dictIO);

  #endif


  /* // Create IOobject for the dictionary */
  /* IOobject dictIO */
  /*   ( */
  /*     dictName, */
  /*     runTime.system(), */
  /*     mesh, */
  /*     IOobject::MUST_READ */
  /*   ); */

  // Read the dictionary from the specified path
  /* dictIO.path() = dictPath; */


  #ifdef OF12
  int iterations = dict.lookupOrDefault<int>("iters", 10);
  int boundaryNormalFreq = dict.lookupOrDefault<int>("boundaryNormalFreq", 1);
  scalar smoothingFactor = dict.lookupOrDefault<scalar>("smoothFactor", 1.0);
  scalar preserveBoundaryLayer = dict.lookupOrDefault<scalar>("preserveBoundaryLayer", 0.0);
  wordReList boundaryNormalPatches = dict.lookup("boundaryNormalPatches");
  #else
  int iterations = dict.getOrDefault<int>("iters", 10);
  int boundaryNormalFreq = dict.getOrDefault<int>("boundaryNormalFreq", 1);
  scalar smoothingFactor = dict.getOrDefault<scalar>("smoothFactor", 1.0);
  scalar preserveBoundaryLayer = dict.getOrDefault<scalar>("preserveBoundaryLayer", 0.0);
  wordReList boundaryNormalPatches = dict.lookup("boundaryNormalPatches");
  #endif
  Info << "iters: " << iterations << endl;
  Info << "smoothFactor: " << smoothingFactor << endl;
  Info << "preserveBoundaryLayer: " << preserveBoundaryLayer << endl;
  Info << "boundaryNormalPatches: " << boundaryNormalPatches << endl;
  if (preserveBoundaryLayer > 0.0)
  {
    WarningInFunction << "Preserve boundary layer functionality is experimental and may not work as expected." << endl;
  }

  // Read set construct info from dictionary
  List<dictionary> constrainedPointsDict(dict.lookup("constrainedPoints"));

  List<word> constrainedNames;
  List<word> constrainedPatchNames;
  List<pointSet*> constrainedPointSets;
  List<word> constraintTypes;
  List<scalar> constraintValues;
  List<vector> constraintDirections;

  for (const dictionary& constrainEntry : constrainedPointsDict)
  {
    const dictionary& constrainedDict = constrainEntry;

    pointSet* ps = nullptr;

    Info << "Checking for type entry in dictionary..." << endl;

    // Read the "type" entry from the dictionary
    if (!constrainedDict.found("type"))
    {
      FatalErrorInFunction << "Missing 'type' entry in dictionary" << endl;
    }
    else
  {
      #ifdef OF12
      word type = constrainedDict.lookup("type");
      #else
      word type = constrainedDict.get<word>("type");
      #endif
      Info << " type: " << type << endl;
      if (type == "patch")
      {
        #ifdef OF12
        word patchName = constrainedDict.lookup("patch");
        #else
        word patchName = constrainedDict.get<word>("patch");
        #endif
        Info << " patch: " << patchName << endl;
        constrainedNames.append(patchName);
        /* label patchID = mesh.boundaryMesh().findPatchID(patchName); */
        /* polyPatch& currentPatch = mesh.boundaryMesh()[patchID]; */
        ps = new pointSet(mesh, patchName, 0);
       forAll(mesh.boundaryMesh()[patchName], faceI)
       {
          face currentFace = mesh.boundaryMesh()[patchName][faceI];
          forAll(currentFace, pointI)
          {
            ps->insert(currentFace[pointI]);
          }
        }
        Info << "Constraining points in patch " << patchName << " with " << ps->size() << " poinst." << endl;
      }
      else {
        #ifdef OF12
        word pointSetName = constrainedDict.lookup("set");
        #else
        word pointSetName = constrainedDict.get<word>("set");
        #endif
        constrainedNames.append(pointSetName);
        ps = new pointSet(mesh, pointSetName);
        Info << "Constraining points in set " << pointSetName << " with " << ps->size() << " poinst." << endl;
      }
    }


    constrainedPointSets.append(ps);
    Info << " number of points: " << constrainedPointSets.last()->size() << endl;

    #ifdef OF12
    word constraintType = constrainedDict.lookup("constraintType");
    #else
    word constraintType = constrainedDict.get<word>("constraintType");
    #endif
    Info << " constraintType: " << constraintType << endl;
    constraintTypes.append(constraintType);

    #ifdef OF12
    scalar constraintValue = constrainedDict.lookupOrDefault<scalar>("value", 1.0);
    #else
    scalar constraintValue = constrainedDict.getOrDefault<scalar>("value", 1.0);
    #endif
    constraintValues.append(constraintValue);

    vector constraintDirection(0, 0, 0);
    if (constraintType == "constDir")
    {
      #ifdef OF12
      constraintDirection = constrainedDict.lookup("direction");
      #else
      constraintDirection = constrainedDict.get<vector>("direction");
      #endif
      Info << " constraintDirection: " << constraintDirection << endl;
    }
    constraintDirections.append(constraintDirection);

    Info << " constraintValue: " << constraintValue << endl;
  }


  const word oldInstance = mesh.pointsInstance();
  #ifdef OF12
  const bool overwrite = args.options().found("overwrite");
  #else
  const bool overwrite = args.found("overwrite");
  #endif

  const polyBoundaryMesh& patches = mesh.boundaryMesh();
  label nPoints = mesh.nPoints();

  // Count how many patches each point belongs to
  labelList patchCount(nPoints, 0);

  forAll(patches, patchI)
  {
      const polyPatch& pp = patches[patchI];
      const labelList& patchPts = pp.meshPoints();

      for (label pt : patchPts)
      {
          patchCount[pt]++;
      }
  }
  
  // labelList boundaryNormalPatchIDs;
  // forAll(boundaryNormalPatches, i)
  // {
  //     word patchName = boundaryNormalPatches[i];
  //     #ifdef OF12
  //     label patchID = mesh.boundaryMesh().findIndex(patchName);
  //     #else
  //     label patchID = mesh.boundaryMesh().findPatchID(patchName);
  //     #endif
  //     if (patchID == -1)
  //     {
  //         FatalErrorInFunction << "Patch " << patchName << " not found in mesh." << endl;
  //     }
  //     boundaryNormalPatchIDs.append(patchID);
  // }
  //
  // Info << "Boundary normal patch IDs: " << boundaryNormalPatchIDs << endl;

  labelHashSet boundaryNormalPatchesSet(mesh.boundaryMesh().patchSet(boundaryNormalPatches));
  labelList boundaryNormalPatchesIDs = boundaryNormalPatchesSet.toc();

  Info << "Boundary normal patch IDs: " << boundaryNormalPatchesIDs << endl;


  // Collect points belonging to multiple patches
  labelList pointsInMultiplePatches;

  forAll(patchCount, ptI)
  {
      if (patchCount[ptI] > 1)
      {
          pointsInMultiplePatches.append(ptI);
          // Info << "Point " << ptI << " appears in "
          //      << patchCount[ptI] << " patches. It will not be moved.\n";
      }
  }

  // // Before starting iterations, we want to find the points that are present in two or more patches
  // // These points will not be moved at all
  // labelList pointsInMultiplePatches;
  //
  // Info << "Checking for points present in multiple patches..." << endl;
  // forAll(mesh.points(), pointI)
  // {
  //   label numberOfPatchesContainingPoint = 0;
  //   forAll(mesh.boundaryMesh(), patchI)
  //   {
  //     polyPatch pp = mesh.boundaryMesh()[patchI];
  //     // Check if the point is in a patch
  //     if (pp.whichPoint(pointI) != -1)
  //     {
  //       numberOfPatchesContainingPoint++;
  //     }
  //   }
  //   if (numberOfPatchesContainingPoint > 1)
  //   {
  //     pointsInMultiplePatches.append(pointI);
  //     Info << "Point " << pointI << " is present in " << numberOfPatchesContainingPoint << " patches. It will not be moved." << endl;
  //   }
  // }

  fvMesh tempFvMesh (mesh);
  

  volScalarField yCell = wallDist::New(tempFvMesh).y();
  volVectorField nCell = wallDist::New(tempFvMesh).n();

  // Check size of n and y, and compare against number of points

  int ySize = yCell.internalField().size();
  int pSize = mesh.points().size();
  Info << "Size of yCell.internalField(): " << ySize << endl;
  Info << "Number of mesh points: " << pSize << endl;
  if (ySize != pSize)
  {
    Info << "Size of wallDist y field does not match number of mesh points!" << endl;
  }
  // As I thought! y and n are fields for cell center values! Need to interpolate to points!
  // Lets use the interpolate function in Foam::meshTools

  // volPointInterpolation vpi(tempFvMesh);
// volPointInterpolation& vpi = volPointInterpolation::New(tempFvMesh)();
// volPointInterpolation& vpi =
//     const_cast<volPointInterpolation&>(volPointInterpolation::New(tempFvMesh));
  // volPointInterpolation& vpi = volPointInterpolation::New(mesh);




  pointScalarField y = volPointInterpolation::New(tempFvMesh).interpolate(yCell);
  pointVectorField n = volPointInterpolation::New(tempFvMesh).interpolate(nCell);

  // Check size of n and y again
  int yPointSize = y.internalField().size();
  if (yPointSize != pSize)
  {
    Info << "Size of interpolated wallDist y field does not match number of mesh points!" << endl;
    exit(1);
  }

  auto t0 = std::chrono::high_resolution_clock::now();

    PackedBoolList vertOnPatch(mesh.nPoints());
    PackedBoolList vertOnExcludePatch(mesh.nPoints());
    labelList wedgePatches;

    // for (const label patchi : mesh.boundaryMesh().patchID())
    forAll(mesh.boundaryMesh(), patchi)
    {
        // Check if patch type is wedge, symmetry, empty or cyclic
        if (mesh.boundaryMesh()[patchi].type() == "wedge"
            || mesh.boundaryMesh()[patchi].type() == "symmetry"
            || mesh.boundaryMesh()[patchi].type() == "empty"
            || mesh.boundaryMesh()[patchi].type() == "cyclic"
            )
        {
            if (mesh.boundaryMesh()[patchi].type() == "wedge")
            {
              wedgePatches.append(patchi);
            }
            continue;
        }
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        const labelList& meshPoints = pp.meshPoints();

        #ifdef OF12
        if (findIndex(boundaryNormalPatchesIDs, patchi) == -1)
        #else
        if (boundaryNormalPatchesIDs.found(patchi))
        #endif
        {
            Info << "Excluding patch " << mesh.boundaryMesh()[patchi].name() << " from boundary normal adjustment." << endl;
            vertOnExcludePatch.set(meshPoints);
        }

        vertOnPatch.set(meshPoints);


    }



  pointField newPoints(mesh.points().size());
  for (label iter = 0; iter < iterations; ++iter) {
    Info << "Iteration " << iter + 1 << endl;
    int movedPoints = 0;
    label constraintCount = 0;
    double timeFind = 0;
    double timeNeighSearch = 0;
    /* Info << "preservedBoundaryLayer: " << preserveBoundaryLayer << endl; */

    // Create wallDist:
    
    // Need to get the fvMesh to use in wallDist! 
    // But we only have polyMesh here.
    // So we create a temporary fvMesh
    // This is a bit of a hack, but it works for now.


    // Loop over all points in the mesh
    forAll(mesh.points(), pointI) {
      vector sumOfNeighbours(0, 0, 0);
      label neighbourCount = 0;

      t0 = std::chrono::high_resolution_clock::now();
      labelList neighboringPointIndices = findNeighboringPoints(mesh, pointI);
      std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t0;

      timeNeighSearch += elapsed.count();

      forAll(neighboringPointIndices, i)
      {
        sumOfNeighbours += mesh.points()[neighboringPointIndices[i]];
        ++neighbourCount;
      }
      vector averagePosition = sumOfNeighbours / neighbourCount;
      vector movement = smoothingFactor * (averagePosition - mesh.points()[pointI]);
      #ifdef OF12
      if (preserveBoundaryLayer and findIndex(pointsInMultiplePatches, pointI) == -1)
      #else
      if (preserveBoundaryLayer and not pointsInMultiplePatches.found(pointI))
      #endif
      {
        if (y.internalField()[pointI] <= preserveBoundaryLayer and y.internalField()[pointI] > SMALL)
        {
          if (boundaryNormalFreq == 0)
          {
            // With no forced boundary normals, we should remove ALL movement close to the boundary:
            movement *= min(1, y[pointI]/preserveBoundaryLayer);
          }
          else
          {
          movement -= (movement & n.internalField()[pointI]) * n.internalField()[pointI] * min(1, 1 - mag(y[pointI]/preserveBoundaryLayer));
          }
        }
      }
      movedPoints++;
      constraintCount++;
      #ifdef OF12
      if (findIndex(pointsInMultiplePatches, pointI) != -1)
      #else
      t0 = std::chrono::high_resolution_clock::now();
      pointsInMultiplePatches.found(pointI);
      elapsed = std::chrono::high_resolution_clock::now() - t0;
      timeFind += elapsed.count();
      if (pointsInMultiplePatches.found(pointI))
      #endif
      {
        movement = vector(0, 0, 0);
        movedPoints--;
        constraintCount++;
      }
      else 
      {
        bool isAlreadyConstrained = false;
        forAll(constrainedNames, i)
        {
          if (constrainedPointSets[i]->found(pointI))
          {
            constraintCount++;
            if (constraintTypes[i] == "yconst")
            {
              movement.y() = constraintValues[i] - mesh.points()[pointI].y();
              constraintCount++;
              if (isAlreadyConstrained)
              {
                movement = vector(0, 0, 0);
              }
              /* isAlreadyConstrained = true; */
            }
            if (constraintTypes[i] == "sphere")
            {
              point onSphere = projectPointOntoSphere(mesh.points()[pointI] + movement, constraintValues[i]);
              movement = onSphere - mesh.points()[pointI];
              constraintCount++;
              if (isAlreadyConstrained)
              {
                movement = vector(0, 0, 0);
              }
              /* isAlreadyConstrained = true; */
            }

            if (constraintTypes[i] == "constX")
            {
              movement.x() = 0;
                /* isAlreadyConstrained = true; */
            }
            if (constraintTypes[i] == "constY")
            {
              movement.y() = 0;
                isAlreadyConstrained = true;
            }
            if (constraintTypes[i] == "constZ")
            {
              movement.z() = 0;
                /* isAlreadyConstrained = true; */
            }

            if (constraintTypes[i] == "constDir")
            {
              // Constrain movement of a specific direction
              /* Info << "Movement before: " << movement << endl; */
              if (isAlreadyConstrained)
              {
                movement = vector(0, 0, 0);
              }
              else
            {
                // Apply a limiter to make a smooth transition between the constrained and unconstrained movement
                /* scalar y0 = 0.32; */
                /* scalar deltay = 0.035; */
                scalar limiter = 1;
                // Between y0 plus minus deltay, the movement is constrained, another deltay, the movement constraint
                // is linearly reduced to zero
                /* if (constraintDirections[i].y() > 0) */
                /* { */
                /*   if (mesh.points()[pointI].y() > y0 - deltay && mesh.points()[pointI].y() < y0 + deltay) */
                /*   { */
                /*     limiter = 1; */
                /*   } */
                /*   else if (mesh.points()[pointI].y() > y0 + deltay) */
                /*   { */
                /*     limiter = 1 - (mesh.points()[pointI].y() - y0 - deltay) / deltay; */
                /*   } */
                /*   else */
                /*   { */
                /*     limiter = 0; */
                /*   } */
                /* } */

                movement -= limiter*(movement & constraintDirections[i]) * constraintDirections[i];
                /* isAlreadyConstrained = true; */
              }
              /* Info << "Movement after: " << movement << endl; */
            }

            if (constraintTypes[i] == "constRadiusXY")
            {
              if (isAlreadyConstrained)
              {
                movement = vector(0, 0, 0);
              }
              else
            {
                /* vector radialUnitVector = mesh.points()[pointI] / mag(mesh.points()[pointI]); */
                /* vector radialMovement = (movement & radialUnitVector) * radialUnitVector; */
                /* movement -= radialMovement; */
                // New method. Project the point onto a sphere with the given radius, and then calculate the movement
                scalar oldRadiusXY = mag(vector(mesh.points()[pointI].x(), mesh.points()[pointI].y(), 0));
                scalar newRadiusXY = mag(vector(mesh.points()[pointI].x() + movement.x(), mesh.points()[pointI].y() + movement.y(), 0));
                scalar scaleFactor = oldRadiusXY / newRadiusXY;
                vector projectedPoint = vector(mesh.points()[pointI].x() + movement.x(), mesh.points()[pointI].y() + movement.y(), 0) * scaleFactor;
                movement.x() = projectedPoint.x() - mesh.points()[pointI].x();
                movement.y() = projectedPoint.y() - mesh.points()[pointI].y();

                
                /* isAlreadyConstrained = true; */

              }
            }

            if (constraintTypes[i] == "fixed")
            {
              movement = vector(0, 0, 0);
              isAlreadyConstrained = true;
              movedPoints--;
            }
          }
        }
        movedPoints++;
      }


      newPoints[pointI] = mesh.points()[pointI] + movement;

    }

    // Update the mesh with the smoothed points
    Info << "Moved " << movedPoints << " points." << endl;
    Info << "Constrained " << constraintCount << " points." << endl;
    mesh.movePoints(newPoints);

    Info << "Time spent in neighbor search: "
       << timeNeighSearch
       << " seconds." << endl;
    Info << "Time spent in find search: "
       << timeFind
       << " seconds." << endl;
    //
    // Adjust the points on the patches to be orthogonal to the internal points
    if ( boundaryNormalFreq > 0 and (iter + 1) % boundaryNormalFreq == 0 )
    {
        Info << "\nAdjusting boundary points to be orthogonal to internal points..." << endl;
        // Info << "Adjusting boundary points to be orthogonal to internal points." << endl;
      forAllConstIter(labelHashSet, boundaryNormalPatchesSet, iter)
      {
          label patchi = iter.key();
          // Info << "Adjusting points on patch " << mesh.boundaryMesh()[patchi].name() << endl;
          const polyPatch& pp = mesh.boundaryMesh()[patchi];
          const labelList& meshPoints = pp.meshPoints();

          for (const label meshPointi : meshPoints)
          {

              // Skip points on excluded patches
              if (vertOnExcludePatch[meshPointi])
              {
                 // Info << "Skipping point " << meshPointi << " on excluded patch" << endl;
                 continue;
               }

              const point& boundaryPoint = mesh.points()[meshPointi];
              const label localPointI = pp.whichPoint(meshPointi);
              // Info << "  Adjusting point " << meshPointi << ": " << boundaryPoint << endl;
              
              // Get the list of points connected to the boundary point
              const labelList& pPoints = mesh.pointPoints()[meshPointi];

              // Find the point that is not on the boundary if there is one
              label internalPointi = -1;
              for (const label pPointi : pPoints)
              {
                  // Info << "    Checking point " << pPointi << endl;
                  if (!vertOnPatch[pPointi])
                  // Info << "       whichPoint: " << pp.whichPoint(pPointi) << endl;
                  // if (pp.whichPoint(pPointi) == -1)
                  {
                      // Info << "    Found internal point " << pPointi << endl;
                      internalPointi = pPointi;
                      break;
                  }
              }
              
              if (internalPointi < 0)
              {
                      // Info << "No internal point found for boundary point " << meshPointi << endl;    
                      // Info << "   coords: " << boundaryPoint << endl;
                      continue;
              }
              else
              {

                  // Calculate the patch normal at the boundary point, by taking the average of the normals of the faces connected to the point
                  vector boundaryNormal = vector::zero;
                  const labelList& faces = pp.pointFaces()[localPointI];
                  /* const auto faces = pp.pointFaces(); */
                  /* const auto faces = pp.localFaces(); */
                  // Info << "    Faces: " << faces << endl;
                  for (const label facei : faces)
                  {
                      /* Info << "    Face " << facei << ": " << faces[facei] << endl; */
                      const vector& faceNormal = pp.faceNormals()[facei];
                      boundaryNormal += faceNormal;
                  }
                  // boundaryNormal.normalise();
                  boundaryNormal /= mag(boundaryNormal);


                  // For now, assumes that the desired normal is (0,-1,0)
                  /* vector boundaryNormal = vector(0, -1, 0); */
                  vector internalPoint = mesh.points()[internalPointi];
                  vector internalToBoundary = boundaryPoint - internalPoint;
                  vector movement = internalToBoundary - dot(internalToBoundary, boundaryNormal) * boundaryNormal;
                  /* point newPoint = internalPoint + internalToBoundary - dot(internalToBoundary, boundaryNormal) * boundaryNormal; */

                  /* mesh.points()[meshPointi] = newPoint; */
                  newPoints[meshPointi] = mesh.points()[meshPointi] - movement;
              }

              
          }
      }
    }
  }
  forAll(wedgePatches, wpi)
  {
    label patchi = wedgePatches[wpi];
    Info << "Constraining wedge patch " << patchi << endl;
    // forAll(newPoints, npointi)
    // {
      // Check if npointi is in the wedge patch:
      // if (mesh.boundaryMesh()[patchi].meshPoints().found(npointi))
      // if (findIndex(mesh.boundaryMesh()[patchi].meshPoints(), npointi))
      const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>(mesh.boundaryMesh()[patchi]);
      const vector n_ = wpp.n();          // wedge normal direction
      forAll(wpp.meshPoints(), wpi)
      {
        label pointi = wpp.meshPoints()[wpi];
        // const vector a_ = wpp.axis();       // wedge axis
        // scalar cosAngle = wpp.cosAngle();   // cos of wedge angle

        // Project new point in the n direction such that the angle with the axis is the same as the wedge:
        vector pointRel = newPoints[pointi] - mesh.points()[pointi];

        // Remove the component in the normal direction

        vector pointRel_n = (pointRel & n_) * n_;

        // Update the new point:
        newPoints[pointi] -= pointRel_n;

      }
    }
    mesh.movePoints(newPoints);

  if (!overwrite) {
    ++runTime;
  } else {
    mesh.setInstance(oldInstance);
  }

  // Write resulting smoothed mesh
  if (Pstream::master()) {
    #ifdef OF12
    Info << "Writing smoothed mesh to time " << runTime.name() << endl;
    #else
    Info << "Writing smoothed mesh to time " << runTime.timeName() << endl;
    #endif
    mesh.write();
  }

  Info << "\nEnd\n" << endl;

  return 0;
}

// ************************************************************************* //
