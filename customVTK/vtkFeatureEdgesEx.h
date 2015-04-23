/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFeatureEdges.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkFeatureEdges - extract boundary, non-manifold, and/or sharp edges from polygonal data
// .SECTION Description
// vtkFeatureEdges is a filter to extract special types of edges from
// input polygonal data. These edges are either 1) boundary (used by
// one polygon) or a line cell; 2) non-manifold (used by three or more
// polygons); 3) feature edges (edges used by two triangles and whose
// dihedral angle > FeatureAngle); or 4) manifold edges (edges used by
// exactly two polygons). These edges may be extracted in any
// combination. Edges may also be "colored" (i.e., scalar values assigned)
// based on edge type. The cell coloring is assigned to the cell data of
// the extracted edges.

// .SECTION Caveats
// To see the coloring of the liens you may have to set the ScalarMode
// instance variable of the mapper to SetScalarModeToUseCellData(). (This
// is only a problem if there are point data scalars.)

// .SECTION See Also
// vtkExtractEdges

#ifndef __vtkFeatureEdgesEx_h
#define __vtkFeatureEdgesEx_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkFeatureEdges.h"

class VTKFILTERSCORE_EXPORT vtkFeatureEdgesEx : public vtkFeatureEdges
{
public:
  vtkTypeMacro(vtkFeatureEdgesEx,vtkFeatureEdges);
  // Description:
  // Construct object with feature angle = 30; all types of edges extracted
  // and colored.
  static vtkFeatureEdgesEx *New();
  vtkIdType GetOldIdFromCurrentID(vtkIdType currentID);
protected:
  vtkFeatureEdgesEx();
  ~vtkFeatureEdgesEx();
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  

  vtkIdList *oldIdList;
private:
  vtkFeatureEdgesEx(const vtkFeatureEdgesEx&);  // Not implemented.
  void operator=(const vtkFeatureEdgesEx&);  // Not implemented.
};

#endif


