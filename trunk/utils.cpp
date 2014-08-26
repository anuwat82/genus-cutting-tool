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

vtkIdType GetCellID(vtkPolyData *polydata, vtkIdType v1,vtkIdType v2,vtkIdType v3)
{
	unsigned short int n1;
	int i, j, tVerts[3];
	vtkIdType *cells, *tVerts2, n2;

	tVerts[0] = v1;
	tVerts[1] = v2;
	tVerts[2] = v3;

	for (i=0; i<3; i++)
	{
		polydata->GetPointCells(tVerts[i], n1, cells);
		for (j=0; j<n1; j++)
		{
			polydata->GetCellPoints(cells[j], n2, tVerts2);
			if ( (tVerts[0] == tVerts2[0] || tVerts[0] == tVerts2[1] ||
				tVerts[0] == tVerts2[2]) &&
				(tVerts[1] == tVerts2[0] || tVerts[1] == tVerts2[1] ||
				tVerts[1] == tVerts2[2]) &&
				(tVerts[2] == tVerts2[0] || tVerts[2] == tVerts2[1] ||
				tVerts[2] == tVerts2[2]) )
			{
				return cells[j];
			}
		}
	}
	return -1;
}


vtkIdType GetCellID(vtkPolyData *polydata, OmMesh &mesh , OmMesh::FaceHandle fh)
{
	OmMesh::HalfedgeHandle _heh[3];
	_heh[0] = mesh.halfedge_handle(fh);
	_heh[1] = mesh.next_halfedge_handle(_heh[0]);
	_heh[2] = mesh.next_halfedge_handle(_heh[1]);
	if (_heh[0].is_valid() && _heh[1].is_valid() && _heh[2].is_valid())
	{
		vtkIdType v1,v2,v3;
		v1 = mesh.to_vertex_handle(_heh[0]).idx();	
		v2 = mesh.to_vertex_handle(_heh[1]).idx();
		v3 = mesh.to_vertex_handle(_heh[2]).idx();
		return GetCellID(polydata,v1,v2,v3);
	}
	else
		return -1;
}