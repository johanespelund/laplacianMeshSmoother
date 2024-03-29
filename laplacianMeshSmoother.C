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
label getClosestBoundaryPoint(const polyMesh& mesh, label pointI) {
  const vectorField& points = mesh.points();
  const point& targetPoint = points[pointI];

  scalar minDistSqr = GREAT; // Initialize with a very large number
  label closestPointIndex = -1;

  // Iterate over all boundary patches
  forAll(mesh.boundaryMesh(), patchI) {
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
      }
    }
  }

  // Return the square root of the shortest squared distance found
  return closestPointIndex;
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

  // Create IOobject for the dictionary
  IOobject dictIO
    (
      dictName,
      runTime.system(),
      mesh,
      IOobject::MUST_READ
    );

  // Read the dictionary from the specified path
  dictIO.path() = dictPath;

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

  List<word> constrainedPointSetNames;
  List<pointSet*> constrainedPointSets;
  List<word> constraintTypes;
  List<scalar> constraintValues;

  for (const namedDictionary& constrainEntry : constrainedPointsDict)
  {
    const dictionary& constrainedDict = constrainEntry.dict();
    word pointSetName = constrainedDict.get<word>("set");
    Info << "Constraining points in " << pointSetName << endl;
    constrainedPointSetNames.append(pointSetName);

    pointSet* ps = new pointSet(mesh, pointSetName);
    constrainedPointSets.append(ps);
    Info << " number of points: " << constrainedPointSets.last()->size() << endl;

    word constraintType = constrainedDict.get<word>("constraintType");
    Info << " constraintType: " << constraintType << endl;
    constraintTypes.append(constraintType);
    scalar constraintValue = constrainedDict.getOrDefault<scalar>("value", 1.0);
    constraintValues.append(constraintValue);
    Info << " constraintValue: " << constraintValue << endl;
  }

  const word oldInstance = mesh.pointsInstance();
  const bool overwrite = args.found("overwrite");

  for (label iter = 0; iter < iterations; ++iter) {
    Info << "Iteration " << iter + 1 << endl;
    pointField newPoints(mesh.points().size());
    int movedPoints = 0;
    label constraintCount = 0;

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

      if (preserveBoundaryLayer > 0.0)
      {
        label closestBoundaryPoint = getClosestBoundaryPoint(mesh, pointI);
        vector distance = mesh.points()[closestBoundaryPoint] - mesh.points()[pointI];
        if ((mag(distance) < preserveBoundaryLayer) )
        {
          vector boundaryNormal = getBoundaryNormal(mesh, closestBoundaryPoint);
          vector parallelComponent = (movement & boundaryNormal) * boundaryNormal;
          movement-= parallelComponent;
          constraintCount++;
        }
      }

      forAll(constrainedPointSetNames, i)
      {
        if (constrainedPointSets[i]->found(pointI))
        {
          if (constraintTypes[i] == "yconst")
          {
            movement.y() = constraintValues[i] - mesh.points()[pointI].y();
            constraintCount++;
          }
          if (constraintTypes[i] == "sphere")
          {
            point onSphere = projectPointOntoSphere(mesh.points()[pointI] + movement, constraintValues[i]);
            movement = onSphere - mesh.points()[pointI];
            constraintCount++;
          }
          if (constraintTypes[i] == "fixed")
          {
            movement = vector(0, 0, 0);
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
