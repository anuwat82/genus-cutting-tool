#include "GIMmodTruncate.h"
#include "utils.h"

GIMmodTruncate::GIMmodTruncate(void)
{
}


GIMmodTruncate::~GIMmodTruncate(void)
{
}
vtkSmartPointer<vtkMutableUndirectedGraph> GIMmodTruncate::GetGraph()
{
	return graph;
}
vtkSmartPointer<vtkMutableUndirectedGraph> GIMmodTruncate::Init( vtkSmartPointer<vtkPolyData> _polydata, 
																 vtkSmartPointer<vtkMutableUndirectedGraph> _collision_graph,
																 geodesic::GeodesicAlgorithmExact *_geo, int geoSourceVertexID)
{
	seedRemoved = false;
	allFaceRemoved =false;
	firstTruncateDone = false;
	branchRemoved = false;
	largestGraphDone = false;
	mesh.clear();
	vtkPolydata2OpenMesh(_polydata,&mesh );	
	this->geodesicExact = _geo;	
	this->polydata = _polydata;	
	this->geodesicSourceVertexID = geoSourceVertexID;
	vtkSmartPointer<vtkMutableUndirectedGraph> _graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	for(vtkIdType i = 0; i < polydata->GetNumberOfPoints(); i++)
    {
		_graph->AddVertex();
    }
  
	vtkSmartPointer<vtkPoints> points = 
    vtkSmartPointer<vtkPoints>::New();
	points->ShallowCopy(polydata->GetPoints());
	_graph->SetPoints(points);

	// Setup the colors array
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("EdgeColors");
	// Setup two colors 
	unsigned char red[3] = {255, 0, 0};
	unsigned char yellow[3] = {255, 255, 0};
	
	// Setup the collision tag array
	vtkSmartPointer<vtkUnsignedCharArray> collisionTag = vtkSmartPointer<vtkUnsignedCharArray>::New();
	collisionTag->SetNumberOfComponents(1);
	collisionTag->SetName("CollisionTag");

	for (OmMesh::EdgeIter e_it=mesh.edges_begin(); e_it!=mesh.edges_end(); ++e_it) 
	{
		OmMesh::EdgeHandle eh = *e_it;
		
		
		OmMesh::HalfedgeHandle heh = mesh.halfedge_handle(eh,0);
		if (!heh.is_valid())
			heh = mesh.halfedge_handle(eh,1);
		OmMesh::HalfedgeHandle heh_prev = mesh.prev_halfedge_handle(heh);
		int a = mesh.to_vertex_handle(heh).idx();
		int b = mesh.to_vertex_handle(heh_prev).idx();
		_graph->AddEdge(a,b);		
		

		colors->InsertNextTupleValue(yellow);
		collisionTag->InsertNextValue(0);

	}
	_graph->GetEdgeData()->AddArray(colors);
	_graph->GetEdgeData()->AddArray(collisionTag);
	this->graph = _graph;
	
	TagSurroundGraphEdges(_collision_graph,1);


	
	return _graph;
}

void GIMmodTruncate::Process()
{
	
	RemoveSeed(false);	
	seedRemoved = true;
	
	RemoveEdgesThatAdjacentOnlyOneFace(false);
	allFaceRemoved = true;
	
	TruncateGraph(false);
	firstTruncateDone = true;

	EliminateBranchPath(false);
	branchRemoved = true;
		
	EliminateUnusedPath(false);
	largestGraphDone = true;

	graph->Modified();
	mesh.garbage_collection();
	polydata->RemoveDeletedCells();
	polydata->BuildLinks();
	
}

void GIMmodTruncate::Step()
{
	if (!seedRemoved)
	{
		//RemoveBoundarySeed(true);
		RemoveSeed(true);
		seedRemoved = true;
	}
	else if (!allFaceRemoved)
	{
	    RemoveEdgesThatAdjacentOnlyOneFace(true);
		
		allFaceRemoved = true;
	}
	else if (!firstTruncateDone)
	{
		if (candidate_nonTagEdges.size() > 0 || candidate_TagEdges.size() > 0)
			throw;
		
		TruncateGraph(true);
		firstTruncateDone = true;
	}
	else if (!branchRemoved)
	{
		
		int numEdgeBefore = 0; graph->GetNumberOfEdges();
		do
		{
			numEdgeBefore = 0; graph->GetNumberOfEdges();
			int numEdge = 0;
			/*
			do
			{
				numEdge = graph->GetNumberOfEdges();
				EliminateSelfCycle();
			}
			while(numEdge > graph->GetNumberOfEdges());
			TruncateGraph();
			*/

			EliminateBranchPath(true);
			
		}
		while (numEdgeBefore > graph->GetNumberOfEdges());
		
		
		branchRemoved = true;
	}
	else if (!largestGraphDone)
	{
		EliminateUnusedPath(true);
		largestGraphDone = true;

	}


}


void GIMmodTruncate::RemoveSeed(bool step)
{
	vtkSmartPointer<vtkUnsignedCharArray> collisionTag = vtkUnsignedCharArray::SafeDownCast( graph->GetEdgeData()->GetArray("CollisionTag"));
	//find surronding edges and face of geodesicSourceVertexID
	OmMesh::VertexHandle source_vh = mesh.vertex_handle(geodesicSourceVertexID);

	
	for (OmMesh::VertexVertexIter vv_it = mesh.vv_iter(source_vh); vv_it.is_valid(); ++vv_it)
	{
		OmMesh::VertexHandle vh = *vv_it;
		vtkIdType edgeid = graph->GetEdgeId(vh.idx(), source_vh.idx());
		if (edgeid < 0)
			throw;
		graph->RemoveEdge(edgeid);
	}
	std::vector<OmMesh::FaceHandle> deletedFace;
	for (OmMesh::VertexFaceIter vf_it = mesh.vf_iter(source_vh); vf_it.is_valid(); ++vf_it)
	{
		OmMesh::FaceHandle fh = *vf_it;
		OmMesh::HalfedgeHandle heh = mesh.halfedge_handle(fh);
		while ( mesh.from_vertex_handle(heh) == source_vh || mesh.to_vertex_handle(heh) == source_vh)
			heh = mesh.next_halfedge_handle(heh);

		if (!mesh.is_boundary(mesh.edge_handle(heh)))
		{
			vtkIdType endID,startID ;
			endID = mesh.to_vertex_handle(heh).idx();
			startID = mesh.from_vertex_handle(heh).idx();
			geodesic::edge_pointer edgeG = geodesicExact->mesh()->find_edge(startID,endID);
			double best_distance;
			geodesic::SurfacePoint p(edgeG);
			geodesicExact->best_source(p,best_distance);

			hedge_data hedata(endID,startID); //reserve input because of opposite hedge

			vtkIdType edgeid = graph->GetEdgeId(endID,startID);
			if (collisionTag->GetValue(edgeid) == 0)
				candidate_nonTagEdges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
			else
				candidate_TagEdges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
			vtkIdType cellID = GetCellID(polydata,mesh,fh);
			polydata->RemoveCellReference(cellID);
			polydata->DeleteCell(cellID);
			//mesh.delete_face( fh,false);
			deletedFace.push_back(fh);
		}
	}

	for (int f = 0 ; f < deletedFace.size(); f++)
	{
		mesh.delete_face(deletedFace[f],false);
	}
	
	if (step)
	{
		mesh.garbage_collection();
		polydata->RemoveDeletedCells();
		polydata->BuildLinks();
	}
}
#if 0
void GIMmodTruncate::RemoveBoundarySeed(bool step)
{
	//remove all boundary edges and their adjacent faces
	//(emulate original gim truncate)

	std::vector<OmMesh::FaceHandle> removeFaces;
	for (OmMesh::EdgeIter e_it=mesh.edges_begin(); e_it!=mesh.edges_end(); ++e_it) 
	{
		OmMesh::EdgeHandle eh = *e_it;
		
		if (mesh.is_boundary(eh))
		{
			
			OmMesh::HalfedgeHandle checkingHeh = mesh.halfedge_handle(eh,0);			
			if (mesh.is_boundary(checkingHeh))
				checkingHeh = mesh.halfedge_handle(eh,1);
			/*
			vtkIdType endID,startID ;
			endID = mesh.to_vertex_handle(checkingHeh).idx();
			startID = mesh.to_vertex_handle( mesh.prev_halfedge_handle(checkingHeh)).idx();
			geodesic::edge_pointer edgeG = geodesicExact->mesh()->find_edge(startID,endID);
			double best_distance;
			geodesic::SurfacePoint p(edgeG);
			geodesicExact->best_source(p,best_distance);
			hedge_data hedata(startID,endID);
			candidate_edges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
			*/
			
			OmMesh::FaceHandle removeFaceHandle = mesh.face_handle(checkingHeh);
			OmMesh::HalfedgeHandle candidate_heh[2];
			
			candidate_heh[0] = mesh.next_halfedge_handle(checkingHeh);
			candidate_heh[1] = mesh.next_halfedge_handle(candidate_heh[0]);
			for (int i = 0 ; i < 2 ; i++)
			{

				OmMesh::FaceHandle oppositeFace = mesh.opposite_face_handle(candidate_heh[i]);
				vtkIdType endID,startID ;
				endID = mesh.to_vertex_handle(candidate_heh[i]).idx();
				startID = mesh.to_vertex_handle( mesh.prev_halfedge_handle(candidate_heh[i])).idx();
				
				if (!mesh.is_boundary(mesh.edge_handle(candidate_heh[i])))
				{					
					geodesic::edge_pointer edgeG = geodesicExact->mesh()->find_edge(startID,endID);
					double best_distance;
					geodesic::SurfacePoint p(edgeG);
					geodesicExact->best_source(p,best_distance);		
					
					
					//OmMesh::Point opfaceCenter;
					//mesh.calc_face_centroid(oppositeFace,opfaceCenter);
					//double best_distance = (opfaceCenter - faceCenter).length();
					
					hedge_data hedata(endID,startID); //reserve input because of opposite hedge
					candidate_edges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
				}

			}
			vtkIdType endID,startID ;
			endID = mesh.to_vertex_handle(checkingHeh).idx();
			startID = mesh.to_vertex_handle( mesh.prev_halfedge_handle(checkingHeh)).idx();
			vtkIdType edgeID =  graph->GetEdgeId(startID,endID);
			graph->RemoveEdge(edgeID);
					
			removeFaces.push_back(removeFaceHandle);			
		}
	}

	for (int i =0 ; i < removeFaces.size(); i++)
	{
		vtkIdType cellID = GetCellID(polydata,mesh,removeFaces[i]);
		polydata->RemoveCellReference(cellID);
		polydata->DeleteCell(cellID);
		mesh.delete_face(removeFaces[i],false);
	}
	
	if (step)
	{
		mesh.garbage_collection();
		polydata->RemoveDeletedCells();
		polydata->BuildLinks();
	}
	seedRemoved = true;
}
#endif
void GIMmodTruncate::GenerateCandidateEdgesAndDeleteEdgeFace(hedge_data &checkingEdge,vtkSmartPointer<vtkUnsignedCharArray> collisionTag)
{
	OmMesh::HalfedgeHandle checkingHeh = mesh.find_halfedge( mesh.vertex_handle(checkingEdge.startID) , mesh.vertex_handle(checkingEdge.endID));
	if (!checkingHeh.is_valid())
	{
		//this edge (of graph) have adjacent 0 edge
		//remove from map -> will be cutpath
		return;	
	}
	else
	{
		//remove face and edge			
		OmMesh::FaceHandle removeFaceHandle = mesh.face_handle(checkingHeh);
		OmMesh::HalfedgeHandle candidate_heh[2];
		OmMesh::Point faceCenter;
		mesh.calc_face_centroid(removeFaceHandle,faceCenter);
		candidate_heh[0] = mesh.next_halfedge_handle(checkingHeh);
		candidate_heh[1] = mesh.next_halfedge_handle(candidate_heh[0]);
		for (int i = 0 ; i < 2 ; i++)
		{
			OmMesh::FaceHandle oppositeFace = mesh.opposite_face_handle(candidate_heh[i]);
			vtkIdType endID,startID ;
			endID = mesh.to_vertex_handle(candidate_heh[i]).idx();
			startID = mesh.to_vertex_handle( mesh.prev_halfedge_handle(candidate_heh[i])).idx();
			vtkIdType edgeID =  graph->GetEdgeId(endID,startID);
			//if (!mesh.is_boundary(mesh.edge_handle(candidate_heh[i])))
			{		
						
				geodesic::edge_pointer edgeG = geodesicExact->mesh()->find_edge(startID,endID);
				double best_distance;
				geodesic::SurfacePoint p(edgeG);
				geodesicExact->best_source(p,best_distance);		
				/*
				OmMesh::Point opfaceCenter;
				mesh.calc_face_centroid(oppositeFace,opfaceCenter);
				double best_distance = (opfaceCenter - faceCenter).length();
				*/
				hedge_data hedata(endID,startID); //reserve input because of opposite hedge
				vtkIdType edgeid = graph->GetEdgeId(endID,startID);
				if (collisionTag->GetValue(edgeid) == 0)
					candidate_nonTagEdges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
				else
					candidate_TagEdges.insert(std::pair<double,hedge_data>(best_distance,hedata ));					
			}
	
		}
		vtkIdType edgeID =  graph->GetEdgeId(checkingEdge.startID,checkingEdge.endID);
		graph->RemoveEdge(edgeID);
		vtkIdType cellID = GetCellID(polydata,mesh,removeFaceHandle);
		polydata->RemoveCellReference(cellID);
		polydata->DeleteCell(cellID);
		mesh.delete_face( removeFaceHandle,false);
	}
}
void GIMmodTruncate::RemoveEdgesThatAdjacentOnlyOneFace(bool step)
{
	vtkSmartPointer<vtkUnsignedCharArray> collisionTag = vtkUnsignedCharArray::SafeDownCast( graph->GetEdgeData()->GetArray("CollisionTag"));
	while (candidate_TagEdges.size() > 0 || candidate_nonTagEdges.size() > 0)
	{
		while (candidate_nonTagEdges.size() > 0)
		{
		
			hedge_data checkingEdge = candidate_nonTagEdges.begin()->second;
			candidate_nonTagEdges.erase(candidate_nonTagEdges.begin());
			GenerateCandidateEdgesAndDeleteEdgeFace(checkingEdge,collisionTag);
		}
		
		while (candidate_TagEdges.size() > 0)
		{
		
			hedge_data checkingEdge = candidate_TagEdges.begin()->second;
			candidate_TagEdges.erase(candidate_TagEdges.begin());
			GenerateCandidateEdgesAndDeleteEdgeFace(checkingEdge,collisionTag);
		}
			
	} //end while
	if (step)
	{
		graph->Modified();
		mesh.garbage_collection();
		polydata->RemoveDeletedCells();
		polydata->BuildLinks();
	}
}

void GIMmodTruncate::TruncateGraph(bool step)
{
	int numRemoveEdges = 0;
	do
	{
		vtkSmartPointer<vtkIdTypeArray> removeEdgesArray = vtkSmartPointer<vtkIdTypeArray>::New();
		numRemoveEdges = removeEdgesArray->GetNumberOfTuples();
		int numEdges = graph->GetNumberOfEdges();
		vtkSmartPointer<vtkEdgeListIterator> eit =vtkSmartPointer<vtkEdgeListIterator>::New();
		graph->GetEdges(eit);

		while (eit->HasNext())
		{
			vtkEdgeType e = eit->Next();
			vtkIdType vid[2] = {e.Source,e.Target};
			for (int i = 0 ; i < 2 ; i++)
			{
								
				int adjnum = graph->GetDegree(vid[i]);
				if (adjnum == 1)
				{
					removeEdgesArray->InsertNextValue(e.Id);
					break;
				}
			}    
		}

		numRemoveEdges = removeEdgesArray->GetNumberOfTuples();
		if (numRemoveEdges > 0)
			graph->RemoveEdges(removeEdgesArray);
	}
	while(numRemoveEdges > 0);
}

void GIMmodTruncate::EliminateSelfCycle(bool step)
{
	
	for (vtkIdType vid = 0; vid < graph->GetNumberOfVertices(); vid++)
	{
		if (graph->GetDegree(vid) == 3)
		{
			vtkIdType adjvid[3];
			vtkSmartPointer<vtkAdjacentVertexIterator> adjacent =vtkSmartPointer<vtkAdjacentVertexIterator>::New();
			graph->GetAdjacentVertices(vid,adjacent);
				
			adjvid[0] = adjacent->Next();
			adjvid[1] = adjacent->Next();
			adjvid[2] = adjacent->Next();
			vtkIdType endofpathvid[3] ;	
			for (int i = 0 ; i < 3 ;i++)
			{
				vtkIdType nextvid = adjvid[i];
				vtkIdType prevvid = vid;
				bool selfLoopFound = false;
				while (graph->GetDegree(nextvid) == 2)
				{
					if (nextvid == adjvid[(i+1)%3])
					{
						//self loop detected!
						//delete edges (vid , adjvid[i]) and (vid,adjvid[(i+1)%3]) 
						vtkIdType e1 = graph->GetEdgeId(vid , adjvid[i]);
						graph->RemoveEdge(e1);
						vtkIdType e2 = graph->GetEdgeId(vid , adjvid[(i+1)%3]);
						graph->RemoveEdge(e2);
						selfLoopFound = true;
						break;

					}
					else if (nextvid == adjvid[(i+2)%3])
					{
						//self loop detected!
						//delete edges (vid , adjvid[i]) and (vid,adjvid[(i+2)%3])
						vtkIdType e1 = graph->GetEdgeId(vid , adjvid[i]);
						graph->RemoveEdge(e1);
						vtkIdType e2 = graph->GetEdgeId(vid , adjvid[(i+2)%3]);
						graph->RemoveEdge(e2);
						selfLoopFound = true;
						break;
					}
					vtkSmartPointer<vtkAdjacentVertexIterator> subadjacent =vtkSmartPointer<vtkAdjacentVertexIterator>::New();
					graph->GetAdjacentVertices(nextvid,subadjacent);
					vtkIdType nextnextvid = subadjacent->Next();
					if (nextnextvid == prevvid)
						nextnextvid = subadjacent->Next();

					prevvid = nextvid;
					nextvid = nextnextvid;

				}
				
				if (selfLoopFound)
					break;
				else
				{
					endofpathvid[i] = nextvid;

				}
			}
			/*
			for (int i = 0 ; i < 3; i++)
			{
				//check for end of path same id or not
				if (endofpathvid[i] != vid && endofpathvid[(i+1)%3] != vid)
				{
					if (endofpathvid[i] == endofpathvid[(i+1)%3])
					{
						//self cycle detected
						//delete edges (vid , adjvid[i])
						vtkIdType e1 = graph->GetEdgeId(vid , adjvid[i]);
						graph->RemoveEdge(e1);
						break;
					}
				}
			}
			*/
		}
	}
}

void GIMmodTruncate::EliminateBranchPath(bool step)
{
	bool *checkedPoint = new bool [graph->GetNumberOfVertices()];
	memset(checkedPoint,0,sizeof(bool)*graph->GetNumberOfVertices());
	int numEdges = 0;
	do
	{
		numEdges = graph->GetNumberOfEdges();
		vtkSmartPointer<vtkEdgeListIterator> eit =vtkSmartPointer<vtkEdgeListIterator>::New();
		graph->GetEdges(eit);
		
		while (eit->HasNext())
		{
			vtkEdgeType e = eit->Next();
			vtkIdType vid[2] = {e.Source,e.Target};
			if (!checkedPoint[vid[0]] || !checkedPoint[vid[1]])
			{
				checkedPoint[vid[0]] = true;
				checkedPoint[vid[1]] = true; 
				vtkSmartPointer<vtkMutableUndirectedGraph> tmpGraph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
				tmpGraph->DeepCopy(graph);
				vtkIdType eid = tmpGraph->GetEdgeId(vid[0],vid[1]);
				tmpGraph->RemoveEdge(eid);
			
				vtkSmartPointer<vtkBoostBreadthFirstSearch> bfs = vtkSmartPointer<vtkBoostBreadthFirstSearch>::New();
				bfs->SetOriginVertex(vid[0]);
				bfs->SetInputData(tmpGraph);
				bfs->Update();
				int numvertbfs = bfs->GetOutput()->GetNumberOfVertices();
				vtkIntArray* level = vtkIntArray::SafeDownCast(bfs->GetOutput()->GetVertexData()->GetArray("BFS"));
				int result = level->GetValue(vid[1]);
				if (result == level->GetDataTypeValueMax())
				{
					//branch detect
					graph->RemoveEdge(e.Id);
					TruncateGraph(step);
					break;
				}
				
			}
		}
	}
	while(numEdges > graph->GetNumberOfEdges());
	delete [] checkedPoint;
}

void GIMmodTruncate::EliminateUnusedPath(bool step)
{
	//elimainate small connectivity
	vtkSmartPointer<vtkBoostExtractLargestComponent> LconnectedComponents = vtkSmartPointer<vtkBoostExtractLargestComponent>::New();
	LconnectedComponents->SetInputData(graph);
	LconnectedComponents->Update();
	vtkGraph* outputGraph = LconnectedComponents->GetOutput(); 
	graph->ShallowCopy(outputGraph);

}

void GIMmodTruncate::GetNeighborCell(vtkPolyData* source_polydata, vtkIdType vid , int level, std::vector<vtkIdType> &storeIDs)
{
	vtkIdType *cellid;
	unsigned short numCell;
	source_polydata->GetPointCells(vid,numCell,cellid);
			
	for (int j = 0 ; j < numCell; j++)
	{
		storeIDs.push_back(cellid[j]);
		if (level > 1)
		{
			vtkIdType npts(0);
			vtkIdType *pts;
			source_polydata->GetCellPoints(cellid[j], npts,pts);
			for (int k = 0 ; k < npts; k++)
			{
				if (pts[k] != vid)
					GetNeighborCell(source_polydata,pts[k],level-1,storeIDs);
			}
			
		}
		
	}
}

void GIMmodTruncate::TagSurroundGraphEdges(vtkSmartPointer<vtkMutableUndirectedGraph> collision_graph,int level)
{
	if (level < 1)
		throw;
	unsigned char red[3] = {255, 0, 0};
	vtkSmartPointer<vtkUnsignedCharArray> edgecolors = vtkUnsignedCharArray::SafeDownCast( graph->GetEdgeData()->GetArray("EdgeColors"));
	vtkSmartPointer<vtkUnsignedCharArray> collisionTag = vtkUnsignedCharArray::SafeDownCast( graph->GetEdgeData()->GetArray("CollisionTag"));
	
	
	vtkSmartPointer<vtkEdgeListIterator> ceit =vtkSmartPointer<vtkEdgeListIterator>::New();
	collision_graph->GetEdges(ceit);


	bool *cellCheck = new bool[polydata->GetPolys()->GetNumberOfCells()];
	memset(cellCheck,0,sizeof(bool)*polydata->GetPolys()->GetNumberOfCells());

	while (ceit->HasNext())
	{
		vtkEdgeType e = ceit->Next();
		vtkIdType vid[2] = {e.Source,e.Target};
		
		for (int i = 0 ; i < 2 ; i++)
		{
			std::vector<vtkIdType> add_cellIDs;
			GetNeighborCell(polydata,vid[i],level,add_cellIDs);
			for (int j = 0 ; j < add_cellIDs.size(); j++)
			{
				if (!cellCheck[add_cellIDs[j]])
				{
					vtkIdType npts;
					vtkIdType *pts;
					
					polydata->GetCellPoints(add_cellIDs[j],npts,pts);
					if (npts != 3)
						throw;
					for (int k = 0 ; k < 3 ;k++)
					{
						vtkIdType edgeid =  graph->GetEdgeId(pts[k%3],pts[(k+1)%3]);
						if (edgeid < 0) 
							throw;
						collisionTag->SetValue(edgeid,1);
						edgecolors->SetTupleValue(edgeid,red);
					}

					cellCheck[add_cellIDs[j]] = true;
				}
			}
		}		
	}

}