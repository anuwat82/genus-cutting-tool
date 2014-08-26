#pragma once 
#include <string>
#include "VTK_Header.h"
#include "OpenMesh_Header.h"
int GetModelFileName(std::string &filename);
void vtkPolydata2OpenMesh(vtkPolyData *polydata, OmMesh *mesh);
vtkIdType GetCellID(vtkPolyData *polydata, vtkIdType v1,vtkIdType v2,vtkIdType v3);
vtkIdType GetCellID(vtkPolyData *polydata, OmMesh &mesh , OmMesh::FaceHandle fh);