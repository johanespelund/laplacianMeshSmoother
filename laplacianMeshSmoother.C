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

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "edgeCollapser.H"
#include "meshTools.H"
#include "Pair.H"
#include "globalIndex.H"
#include "topoSet.H"
#include "processorMeshes.H"
#include "IOdictionary.H"
#include "namedDictionary.H"
#include "pointSet.H"
#include "faceList.H"

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
      if (currentFace.found(pointIndex))
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
  neighboringPointIndices = Foam::uniqueSort(neighboringPointIndices);

  return neighboringPointIndices;
}

int main(int argc, char *argv[]) {

  argList::addNote("Boundary smoothing utility.");
  #include "addOverwriteOption.H"
  argList::addOption("dict", "file", "Alternative refineMeshDict");
  argList::noFunctionObjects(); // Never use function objects

  #include "setRootCase.H"
  #include "createTime.H"
  #include "createPolyMesh.H"

  // Obtain dictPath here for messages
  fileName dictPath = args.getOrDefault<fileName>("dict", "");

  // Dictionary to control refinement
  dictionary dict;
  const word dictName("laplacianSmoothDict");
  #include "setSystemMeshDictionaryIO.H"

  Info << "Reading dictionary " << dictName << " from " << dictPath << endl;

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

  // Read the dictionary
  dict = IOdictionary(dictIO);

  int iterations = dict.getOrDefault<int>("iters", 10);
  Info << "iters: " << iterations << endl;

  scalar smoothingFactor = dict.getOrDefault<scalar>("smoothFactor", 1.0);
  Info << "smoothFactor: " << smoothingFactor << endl;

  scalar preserveBoundaryLayer = dict.getOrDefault<scalar>("preserveBoundaryLayer", 0.0);
  Info << "preserveBoundaryLayer: " << preserveBoundaryLayer << endl;
  if (preserveBoundaryLayer > 0.0)
  {
    WarningInFunction << "Preserve boundary layer functionality is experimental and may not work as expected." << endl;
  }

  // Read set construct info from dictionary
  List<namedDictionary> constrainedPointsDict(dict.lookup("constrainedPoints"));

  List<word> constrainedNames;
  List<word> constrainedPatchNames;
  List<pointSet*> constrainedPointSets;
  List<word> constraintTypes;
  List<scalar> constraintValues;
  List<vector> constraintDirections;

  for (const namedDictionary& constrainEntry : constrainedPointsDict)
  {
    const dictionary& constrainedDict = constrainEntry.dict();

    pointSet* ps = nullptr;

    Info << "Checking for type entry in dictionary..." << endl;

    // Read the "type" entry from the dictionary
    if (!constrainedDict.found("type"))
    {
      FatalErrorInFunction << "Missing 'type' entry in dictionary" << endl;
    }
    else
  {
      word type = constrainedDict.get<word>("type");
      Info << " type: " << type << endl;
      if (type == "patch")
      {
        word patchName = constrainedDict.get<word>("patch");
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
        word pointSetName = constrainedDict.get<word>("set");
        constrainedNames.append(pointSetName);
        ps = new pointSet(mesh, pointSetName);
        Info << "Constraining points in set " << pointSetName << " with " << ps->size() << " poinst." << endl;
      }
    }


    constrainedPointSets.append(ps);
    Info << " number of points: " << constrainedPointSets.last()->size() << endl;

    word constraintType = constrainedDict.get<word>("constraintType");
    Info << " constraintType: " << constraintType << endl;
    constraintTypes.append(constraintType);

    scalar constraintValue = constrainedDict.getOrDefault<scalar>("value", 1.0);
    constraintValues.append(constraintValue);

    vector constraintDirection(0, 0, 0);
    if (constraintType == "constDir")
    {
      constraintDirection = constrainedDict.get<vector>("direction");
      Info << " constraintDirection: " << constraintDirection << endl;
    }
    constraintDirections.append(constraintDirection);

    Info << " constraintValue: " << constraintValue << endl;
  }

  const word oldInstance = mesh.pointsInstance();
  const bool overwrite = args.found("overwrite");

  for (label iter = 0; iter < iterations; ++iter) {
    Info << "Iteration " << iter + 1 << endl;
    pointField newPoints(mesh.points().size());
    int movedPoints = 0;
    label constraintCount = 0;
    /* Info << "preservedBoundaryLayer: " << preserveBoundaryLayer << endl; */

    // Loop over all points in the mesh
    forAll(mesh.points(), pointI) {
      vector sumOfNeighbours(0, 0, 0);
      label neighbourCount = 0;

      labelList neighboringPointIndices = findNeighboringPoints(mesh, pointI);

      forAll(neighboringPointIndices, i)
      {
        sumOfNeighbours += mesh.points()[neighboringPointIndices[i]];
        ++neighbourCount;
      }
      vector averagePosition = sumOfNeighbours / neighbourCount;
      vector movement = smoothingFactor * (averagePosition - mesh.points()[pointI]);
      int numberOfClosePatches = 0;
      label closestBoundaryPoint = -1;

      if (preserveBoundaryLayer == 0)
      {
        forAll(mesh.boundaryMesh(), patchI)
        {
          if (preserveBoundaryLayer > 0.0)
          {
            closestBoundaryPoint = getClosestBoundaryPoint(mesh, pointI, preserveBoundaryLayer, patchI);
          } else {
            // Check if the point is in a patch
            closestBoundaryPoint = -1;
            forAll(mesh.boundaryMesh()[patchI], faceI)
            {
              face currentFace = mesh.boundaryMesh()[patchI][faceI];
              if (currentFace.found(pointI))
              {
                /* Info << "Found a point on a boundary patch!" << endl; */
                closestBoundaryPoint = pointI;
                break;
              }
            }
          }

          if (closestBoundaryPoint != -1)
          {
            vector distance = mesh.points()[closestBoundaryPoint] - mesh.points()[pointI];
            /* Info << "Distance to boundary: " << mag(distance) << endl; */
            if (mag(distance) <= preserveBoundaryLayer)
            {
              numberOfClosePatches++;
              vector boundaryNormal = getBoundaryNormal(mesh, closestBoundaryPoint);
              vector parallelComponent = (movement & boundaryNormal) * boundaryNormal;
              movement-= parallelComponent;
              /* Info << "Constrained movement to boundary normal." << endl; */
            }
          }
          /* else */
          /* { */
          /*   movement = vector(0, 0, 0); */
          /*   movedPoints--; */
          /*   constraintCount++; */
          /* } */
        }
        movedPoints++;
        constraintCount++;
        if (numberOfClosePatches > 1){
          /* Info << "Found a point which is close to two patches! Not moving at all." << endl; */
          movement = vector(0, 0, 0);
          movedPoints--;
        }
      }
      /* label closestBoundaryPoint = getClosestBoundaryPoint(mesh, pointI, preserveBoundaryLayer); */
      /* if (closestBoundaryPoint != -1) */
      /* { */
      /*   vector distance = mesh.points()[closestBoundaryPoint] - mesh.points()[pointI]; */
      /*     vector boundaryNormal = getBoundaryNormal(mesh, closestBoundaryPoint); */
      /*     vector parallelComponent = (movement & boundaryNormal) * boundaryNormal; */
      /*     scalar x = mag(distance); */
      /*     scalar y = 8/7 - x/(7*preserveBoundaryLayer);  //x/(3*preserveBoundaryLayer) - 1/3; */
      /*   // TODO: Restrict movement in two directions for points close to the boundary */
      /*     movement-= parallelComponent*min(1,max(0, y)); */
      /*     movedPoints++; */
      /*     constraintCount++; */
      /* } */
      /* else */
      /* { */
      /*   movement = vector(0, 0, 0); */
      /*   movedPoints--; */
      /*   constraintCount++; */
      /* } */

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

      newPoints[pointI] = mesh.points()[pointI] + movement;
      movedPoints++;

    }

    // Update the mesh with the smoothed points
    Info << "Moved " << movedPoints << " points." << endl;
    Info << "Constrained " << constraintCount << " points." << endl;
    mesh.movePoints(newPoints);
  }

  if (!overwrite) {
    ++runTime;
  } else {
    mesh.setInstance(oldInstance);
  }

  // Write resulting smoothed mesh
  if (Pstream::master()) {
    Info << "Writing smoothed mesh to time " << runTime.timeName() << endl;
    mesh.write();
  }

  Info << "\nEnd\n" << endl;

  return 0;
}

// ************************************************************************* //
