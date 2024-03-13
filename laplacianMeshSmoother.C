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

// Function to get point IDs belonging to a patch
List<label> getPatchPointIDs(const polyMesh &mesh, word patch_name) {
  
  label patchID = mesh.boundaryMesh().findPatchID(patch_name);

  List<int> patchPointIDs;

  // Access patch faces
  const polyPatch &patch = mesh.boundaryMesh()[patchID];
  Info << patch_name << " : " << patchID << endl;

  // Iterate over faces of the patch
  forAll(patch, faceID) {
    const face &meshFace = patch[faceID];

    Info << faceID;

    // Access the points of the face
    const List<int> &facePoints = meshFace;

    Info << " " << facePoints << endl;

    // Append point IDs to the list
    patchPointIDs.append(facePoints);
  }

  labelList result;
  Info << "patchPointIDs min/max:" << min(patchPointIDs) << " "
       << max(patchPointIDs) << endl;

  // Remove duplicates
  result = Foam::uniqueSort(patchPointIDs);
  /* Foam::sort(patchPointIDs); */
  Info << "result min/max:" << min(result) << " " << max(result) << endl;
  Info << "result: " << result << endl;

  return result;
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


// A class to holde the constrained point set:
// Input: pointSetName, constraintType, constraintValue
// Members: pointSetName, pointSet, constraintType, constraintValue

/* class constrainedPointSet */
/* { */
/*   public: */
/*     word pointSetName; */
/*     pointSet ps; */
/*     word constraintType; */
/*     scalar constraintValue; */

/*     constrainedPointSet(const polyMesh &mesh, word pointSetName, word constraintType, scalar constraintValue) */
/*     { */
/*       this->pointSetName = pointSetName; */
/*       this->ps = pointSet(mesh, pointSetName); */
/*       this->constraintType = constraintType; */
/*       this->constraintValue = constraintValue; */
/*     } */
/* }; */

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

  /* argList::addOption("iters", "int", "Number of smoothing iterations."); */
  /* argList::addOption("smooth-factor", "scalar", "Smoothing factor."); */
  /* /1* argList::addOption("patch", "name", "Name of the patch to smooth."); *1/ */
  /* argList::addOption("points", "name", "Name of the pointSet to smooth."); */
  /* argList::addOption("constant-points", "name", "Name of the pointSet to which should not be moved."); */
  /* argList::addOption("constrain-type", "word", "Constraint type (yconst or rconst)."); */
  /* argList::addOption("value", "scalar", "Value for the constraint (radius or y-coordinate)."); */
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
      Info << " number of points: " << ps->size() << endl;

      word constraintType = constrainedDict.get<word>("constraintType");
      Info << " constraintType: " << constraintType << endl;
      constraintTypes.append(constraintType);
      scalar constraintValue = constrainedDict.getOrDefault<scalar>("value", 1.0);
      constraintValues.append(constraintValue);
      Info << " constraintValue: " << constraintValue << endl;
    }

    /* // Read constrainedPoints */
    /* const dictionary& constrainedPointsDict = dict.subDict("constrainedPoints"); */

    /* Info << "constrainedPoints:" << endl; */
  /* Info << constrainedPointsDict << endl; */
  // Read the dictionary. It looks like this:
/* intenalPoints internalPoints; */
/* iters 10; */
/* smoothFactor 1.0; */

/* constrainedPoints */
/* ( */
/*   { */
/*     set wallSmoothPoints; */
/*     constraintType sphere; */
/*     value $R; */
/*   } */
/*   { */
/*     set bottomSmoothPoints; */
/*     constraintType yconst; */
/*     value #eval {$R - $H_G}; */
/*   } */
/*   { */
/*     set fixedPoints; */
/*     constraintType fixed; */
/*   } */
/* ); */


  const word oldInstance = mesh.pointsInstance();
  const bool overwrite = args.found("overwrite");

  /* label iterations = args.found("iters") */
  /*                        ? readLabel(args["iters"]) */
  /*                        : 5; // Default to 10 iterations if not provided */

  /* scalar smoothingFactor = */
  /*     args.found("smooth-factor") */
  /*         ? readScalar(args["smooth-factor"]) */
  /*         : 1.0; // Default smoothing factor to 0.1 if not provided */

  /* word patch_name = args.found("patch") ? word(args["patch"]) : ""; */
  /* scalar value = */
  /*     args.found("value") */
  /*         ? readScalar(args["value"]) */
  /*         : 1; // Default smoothing factor to 0.1 if not provided */
  /* word contrain_type = args.found("constrain-type") ? word(args["constrain-type"]) : ""; */


  /* Info << "Starting Laplacian smoothing with factor: " << smoothingFactor */
  /*      << " and iterations: " << iterations << endl; */

  /* word pointSet_name = ""; */
  /* word const_pointSet_name = ""; */

  /* if (args.found("points")) { */
  /*   pointSet_name = word(args["points"]); */
  /*   Info << "Reading points " << pointSet_name << endl; */
  /* } else { */
  /*   FatalErrorInFunction << "No pointSet specified." << abort(FatalError); */
  /* } */

  /* if (args.found("constant-points")) */
  /* { */
  /*   const_pointSet_name = word(args["constant-points"]); */
  /*   Info << "Found constant points "  << const_pointSet_name << endl; */
  /* } */
  /* else { */
  /*   // Need to set the name to something, otherwise the pointSet will be empty */
  /*   const_pointSet_name = pointSet_name; */
  /* } */

  /* pointSet const_points(mesh, const_pointSet_name); */ 
  /* /1* Info << const_points << endl; *1/ */


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
        newPoints[pointI] =
          mesh.points()[pointI] +
          smoothingFactor * (averagePosition - mesh.points()[pointI]);
        movedPoints++;

          forAll(constrainedPointSetNames, i)
          {
            if (constrainedPointSets[i]->found(pointI))
            {
              /* Info << "Constraining point " << pointI << " in " << constrainedPointSetNames[i] << endl; */
              if (constraintTypes[i] == "yconst")
              {
                /* Info << "Constraining y-coordinate to " << constraintValues[i] << endl; */
                newPoints[pointI].y() = constraintValues[i];
                constraintCount++;
              }
              else if (constraintTypes[i] == "sphere")
              {
                /* Info << "Projecting point onto sphere with radius " << constraintValues[i] << endl; */
                newPoints[pointI] = projectPointOntoSphere(newPoints[pointI], constraintValues[i]);
                constraintCount++;
              }
              else if (constraintTypes[i] == "fixed")
              {
                /* Info << "Keeping point fixed." << endl; */
                newPoints[pointI] = mesh.points()[pointI];
            movedPoints--;
              }
            }
          }
    }

    Info << "Moved " << movedPoints << " points." << endl;
    /* Info << "Constrained " << constraintCount << " points." << endl; */

    // Update the mesh with the smoothed points
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
