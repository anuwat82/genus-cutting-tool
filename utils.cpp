#include "utils.h"
#include <Windows.h>
#include <string>


std::string GetFileExtension(const std::string& FileName)
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}


int GetFileName(std::string &filename ,LPTSTR dialogTitle ,const char* filter ,bool save)
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
	ofn.lpstrFilter = filter;
	ofn.nFilterIndex = 1;
	ofn.lpstrFileTitle = NULL;
	ofn.nMaxFileTitle = 0;
	ofn.lpstrInitialDir = NULL;
	ofn.lpstrTitle = dialogTitle;
	if (!save)
		ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
	else
		ofn.Flags = OFN_SHOWHELP | OFN_OVERWRITEPROMPT; 
	// Display the Open dialog box. 
	if (!save)
	{
		if (GetOpenFileName(&ofn)==TRUE) 
		{
			filename = std::string(ofn.lpstrFile);
			return 1;
		}		
	}
	else
	{
		if (GetSaveFileName(&ofn) == TRUE)
		{
			filename = std::string(ofn.lpstrFile);
			return 1;
		}
	}
	
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


vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, vtkIdType ptId)
{

	vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New(); //output
  //get all cells that vertex 'seed' is a part of
  vtkSmartPointer<vtkIdList> cellIdList =
      vtkSmartPointer<vtkIdList>::New();
  mesh->GetPointCells(ptId, cellIdList);


  //loop through all the cells that use the seed point
  for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
    {

    vtkCell* cell = mesh->GetCell(cellIdList->GetId(i));
    //cout << "The cell has " << cell->GetNumberOfEdges() << " edges." << endl;

    //if the cell doesn't have any edges, it is a line
    if(cell->GetNumberOfEdges() <= 0)
      {
      continue;
      }

    for(vtkIdType e = 0; e < cell->GetNumberOfEdges(); e++)
      {
      vtkCell* edge = cell->GetEdge(e);

      vtkIdList* pointIdList = edge->GetPointIds();
    
      if(pointIdList->GetId(0) == ptId || pointIdList->GetId(1) == ptId)
        {
        if(pointIdList->GetId(0) == ptId)
          {
          connectedVertices->InsertUniqueId(pointIdList->GetId(1));
		 
          }
        else
          {
          connectedVertices->InsertUniqueId(pointIdList->GetId(0));
          }
        }
      }


    }
  
  return connectedVertices;

} 


void CutMesh(OmMesh &mesh, OmMesh::HalfedgeHandle &he0 ,OmMesh::HalfedgeHandle &he1,OmMesh::VertexHandle newVH)
{
	if (mesh.to_vertex_handle(he0) != mesh.from_vertex_handle(he1))
		throw;
	
	OmMesh::VertexHandle oldVH = mesh.to_vertex_handle(he0);
	OmMesh::HalfedgeHandle oph0 = mesh.opposite_halfedge_handle(he0);
	OmMesh::HalfedgeHandle oph1 = mesh.opposite_halfedge_handle(he1);
	bool boundary0 = mesh.is_boundary(oph0);
	bool boundary1 = mesh.is_boundary(oph1);
	OmMesh::HalfedgeHandle prevh0 = mesh.prev_halfedge_handle(he0);
	OmMesh::HalfedgeHandle prevh1 = mesh.prev_halfedge_handle(he1);
	OmMesh::HalfedgeHandle nh0 = mesh.next_halfedge_handle(he0);
	OmMesh::HalfedgeHandle nh1 = mesh.next_halfedge_handle(he1);

	bool val2 = false;
	if (nh0 == he1)
		val2 = true;

	//add the new edge
	OmMesh::HalfedgeHandle new_e0 = mesh.new_edge(mesh.from_vertex_handle(he0), newVH);
	OmMesh::HalfedgeHandle new_e1 = mesh.new_edge(newVH ,mesh.to_vertex_handle(he1));
	mesh.set_face_handle(new_e0, mesh.face_handle(prevh0));
	mesh.set_face_handle(new_e1, mesh.face_handle(prevh1));
	mesh.set_boundary(mesh.opposite_halfedge_handle(new_e0));
	mesh.set_boundary(mesh.opposite_halfedge_handle(new_e1));
	mesh.set_next_halfedge_handle(prevh0,new_e0);
	if (!val2)
		mesh.set_next_halfedge_handle(prevh1,new_e1);
	
	if (!val2)
		mesh.set_next_halfedge_handle(new_e0,nh0);
	else
		mesh.set_next_halfedge_handle(new_e0,new_e1);
	mesh.set_next_halfedge_handle(new_e1,nh1);

	OmMesh::HalfedgeHandle start_heh = nh0;
	OmMesh::HalfedgeHandle next_heh = start_heh;
	while (next_heh != he1)
	{
		mesh.set_vertex_handle( mesh.opposite_halfedge_handle(next_heh), newVH);
		next_heh = mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(next_heh));
	}
	if (!boundary0)
		mesh.set_boundary(he0);
	if (!boundary1)
		mesh.set_boundary(he1);
	mesh.set_halfedge_handle( newVH, new_e0 );
	mesh.adjust_outgoing_halfedge( newVH );

	mesh.set_halfedge_handle( oldVH, oph1 );
	mesh.adjust_outgoing_halfedge( oldVH );

	if (!boundary0)
		mesh.set_next_halfedge_handle(he0,he1); //boundary
	


	/*
	HalfedgeHandle h0 = halfedge_handle(_eh, 0);
	HalfedgeHandle h1 = halfedge_handle(_eh, 1);
  
	VertexHandle vfrom = from_vertex_handle(h0);

	HalfedgeHandle ph0 = prev_halfedge_handle(h0);
	HalfedgeHandle ph1 = prev_halfedge_handle(h1);
  
	HalfedgeHandle nh0 = next_halfedge_handle(h0);
	HalfedgeHandle nh1 = next_halfedge_handle(h1);
  
	bool boundary0 = is_boundary(h0);
	bool boundary1 = is_boundary(h1);
  
	//add the new edge
	HalfedgeHandle new_e = new_edge(from_vertex_handle(h0), _vh);
  
	//fix the vertex of the opposite halfedge
	set_vertex_handle(h1, _vh);
  
	//fix the halfedge connectivity
	set_next_halfedge_handle(new_e, h0);
	set_next_halfedge_handle(h1, opposite_halfedge_handle(new_e));
  
	set_next_halfedge_handle(ph0, new_e);
	set_next_halfedge_handle(opposite_halfedge_handle(new_e), nh1);
  
	//  set_prev_halfedge_handle(new_e, ph0);
	//  set_prev_halfedge_handle(opposite_halfedge_handle(new_e), h1);
  
	//  set_prev_halfedge_handle(nh1, opposite_halfedge_handle(new_e));
	//  set_prev_halfedge_handle(h0, new_e);
  
	if (!boundary0)
	{
	set_face_handle(new_e, face_handle(h0));
	}
	else
	{
	set_boundary(new_e);
	}
  
	if (!boundary1)
	{
	set_face_handle(opposite_halfedge_handle(new_e), face_handle(h1));
	}
	else
	{
	set_boundary(opposite_halfedge_handle(new_e));
	}

	set_halfedge_handle( _vh, h0 );
	adjust_outgoing_halfedge( _vh );
  
	if (halfedge_handle(vfrom) == h0)
	{
	set_halfedge_handle(vfrom, new_e);
	adjust_outgoing_halfedge( vfrom );
	}
	*/
}