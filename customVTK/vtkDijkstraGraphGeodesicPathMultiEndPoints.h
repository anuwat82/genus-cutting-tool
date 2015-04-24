/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDijkstraGraphGeodesicPath.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkDijkstraGraphGeodesicPath - Dijkstra algorithm to compute the graph geodesic.
// .SECTION Description
// Takes as input a polygonal mesh and performs a single source shortest
// path calculation. Dijkstra's algorithm is used. The implementation is
// similar to the one described in Introduction to Algorithms (Second Edition)
// by Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and
// Cliff Stein, published by MIT Press and McGraw-Hill. Some minor
// enhancement are added though. All vertices are not pushed on the heap
// at start, instead a front set is maintained. The heap is implemented as
// a binary heap. The output of the filter is a set of lines describing
// the shortest path from StartVertex to EndVertex.
//
// .SECTION Caveats
// The input polydata must have only triangle cells.
//
// .SECTION Thanks
// The class was contributed by Rasmus Paulsen.
// www.imm.dtu.dk/~rrp/VTK . Also thanks to Alexandre Gouaillard and Shoaib
// Ghias for bug fixes and enhancements.

#ifndef vtkDijkstraGraphGeodesicPathMultiEndPoints_h
#define vtkDijkstraGraphGeodesicPathMultiEndPoints_h

#include "vtkFiltersModelingModule.h" // For export macro
#include "vtkDijkstraGraphGeodesicPath.h"

class VTKFILTERSMODELING_NO_EXPORT vtkDijkstraGraphGeodesicPathMultiEndPoints :
                           public vtkDijkstraGraphGeodesicPath
{
public:

  // Description:
  // Instantiate the class
  static vtkDijkstraGraphGeodesicPathMultiEndPoints *New();

  // Description:
  // Standard methods for printing and determining type information.
  vtkTypeMacro(vtkDijkstraGraphGeodesicPathMultiEndPoints,vtkDijkstraGraphGeodesicPath);
  

protected:
  vtkDijkstraGraphGeodesicPathMultiEndPoints();
  ~vtkDijkstraGraphGeodesicPathMultiEndPoints(); 

private:
  vtkDijkstraGraphGeodesicPathMultiEndPoints(const vtkDijkstraGraphGeodesicPathMultiEndPoints&);  // Not implemented.
  void operator=(const vtkDijkstraGraphGeodesicPathMultiEndPoints&);  // Not implemented.

};

#endif

