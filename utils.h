#pragma once 
#include <string>
#include "VTK_Header.h"
#include "OpenMesh_Header.h"

int GetFileName(std::string &filename,const char* filter,bool save = false);
std::string GetFileExtension(const std::string& FileName);

void vtkPolydata2OpenMesh(vtkPolyData *polydata, OmMesh *mesh);
vtkIdType GetCellID(vtkPolyData *polydata, vtkIdType v1,vtkIdType v2,vtkIdType v3);
vtkIdType GetCellID(vtkPolyData *polydata, OmMesh &mesh , OmMesh::FaceHandle fh);
//bool	  isBoundaryEdge(vtkPolyData *polydata, vtkIdType v1,vtkIdType v2);
vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, vtkIdType vid);

void CutMesh(OmMesh &mesh, OmMesh::HalfedgeHandle &he0 ,OmMesh::HalfedgeHandle &he1,OmMesh::VertexHandle newVH);