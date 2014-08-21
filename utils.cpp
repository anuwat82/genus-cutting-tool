#include "utils.h"
#include <Windows.h>
#include <string>
int GetModelFileName(std::string &filename)
{
	OPENFILENAME ofn;       // common dialog box structure
	char szFile[260];       // buffer for file name
	HWND hwnd;              // owner window
	HANDLE hf;              // file handle

	// Initialize OPENFILENAME
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = NULL;
	ofn.lpstrFile = szFile;
	// Set lpstrFile[0] to '\0' so that GetOpenFileName does not 
	// use the contents of szFile to initialize itself.
	ofn.lpstrFile[0] = '\0';
	ofn.nMaxFile = sizeof(szFile);
	ofn.lpstrFilter = "All\0*.*\0PLY\0*.ply\0OFF\0*.off\0";
	ofn.nFilterIndex = 1;
	ofn.lpstrFileTitle = NULL;
	ofn.nMaxFileTitle = 0;
	ofn.lpstrInitialDir = NULL;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

	// Display the Open dialog box. 

	if (GetOpenFileName(&ofn)==TRUE) 
	{
		filename = std::string(ofn.lpstrFile);
		return 1;
	}
	else
		return 0;
		
}


void vtkPolydata2OpenMesh(vtkPolyData *polydata, OmMesh *mesh)
{
	vtkIdType numOfPoints =polydata->GetNumberOfPoints();
	vtkIdType numOfFaces = polydata->GetPolys()->GetNumberOfCells();
	OmMesh::VertexHandle *vh = new OmMesh::VertexHandle[numOfPoints];
	for (vtkIdType vid = 0  ; vid < numOfPoints; vid++)
	{
		vh[vid] = mesh->add_vertex(OmMesh::Point(polydata->GetPoint(vid)));
	}

	polydata->GetPolys()->InitTraversal();
	vtkIdType npts;
	vtkIdType *pointID;
	while(polydata->GetPolys()->GetNextCell(npts,pointID) != 0)
	{
		if (npts != 3)
			throw;
		mesh->add_face(vh[*pointID],vh[*(pointID+1)],vh[*(pointID+2)]);
	}

	if (mesh->n_faces() != numOfFaces)
		throw;
	if (mesh->n_vertices() != numOfPoints)
		throw;
}