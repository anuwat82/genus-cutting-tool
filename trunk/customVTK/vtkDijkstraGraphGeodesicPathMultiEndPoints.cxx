/*=========================================================================

  Program:   Visualization Toolkit
  Module:  vtkDijkstraGraphGeodesicPath.cxx
  Language:  C++
  Date:    $Date$
  Version:   $Revision$

  Made by Rasmus Paulsen
  email:  rrp(at)imm.dtu.dk
  web:    www.imm.dtu.dk/~rrp/VTK

  This class is not mature enough to enter the official VTK release.
=========================================================================*/
#include "vtkDijkstraGraphGeodesicPathMultiEndPoints.h"

#include "vtkCellArray.h"
//#include "vtkDijkstraGraphInternals.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"


vtkStandardNewMacro(vtkDijkstraGraphGeodesicPathMultiEndPoints);

// Construct object with feature angle = 30; all types of edges, except
// manifold edges, are extracted and colored.
vtkDijkstraGraphGeodesicPathMultiEndPoints::vtkDijkstraGraphGeodesicPathMultiEndPoints()
{
	
}

vtkDijkstraGraphGeodesicPathMultiEndPoints::~vtkDijkstraGraphGeodesicPathMultiEndPoints()
{
	
}

