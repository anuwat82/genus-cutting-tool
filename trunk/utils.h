#pragma once 
#include <string>
#include "VTK_Header.h"
#include "OpenMesh_Header.h"
int GetModelFileName(std::string &filename);
void vtkPolydata2OpenMesh(vtkPolyData *polydata, OmMesh *mesh);