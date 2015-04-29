#include "GIMmodTruncate.h"
#include "customVTK\vtkFeatureEdgesEx.h"
#include "customVTK\vtkDijkstraGraphGeodesicPathMultiStartEndPoints.h"
#include "utils.h"

GIMmodTruncate::GIMmodTruncate(void)
{
	seedRemoved = false;
	allFaceRemoved =false;
	firstTruncateDone = false;
	shortenRingsDone = false;
	removePathsBetweenBoundariesFromGraph = false;
	timeConsumed = 0;
}


GIMmodTruncate::~GIMmodTruncate(void)
{
}
vtkSmartPointer<vtkMutableUndirectedGraph> GIMmodTruncate::GetGraph()
{
	return graph;
}

bool GIMmodTruncate::isReadyToCut()
{
	return (firstTruncateDone ||shortenRingsDone);
}
vtkSmartPointer<vtkMutableUndirectedGraph> GIMmodTruncate::Init( vtkSmartPointer<vtkPolyData> _polydata, 
																 vtkSmartPointer<vtkMutableUndirectedGraph> _collision_graph,
																 geodesic::GeodesicAlgorithmExact *_geo, int geoSourceVertexID)
{
	timeConsumed = 0;
	clock_t start = clock();
	seedRemoved = false;
	allFaceRemoved =false;
	firstTruncateDone = false;
	shortenRingsDone = false;
	removePathsBetweenBoundariesFromGraph = false;
	mesh.clear();
	vtkPolydata2OpenMesh(_polydata,&mesh );	

	
	this->geodesicExact = _geo;	
	
	this->polydata = _polydata;	
	this->original_polydata = vtkSmartPointer<vtkPolyData>::New();
	this->original_polydata->DeepCopy(_polydata);
	this->original_polydata->BuildLinks();
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
	colors->SetNumberOfComponents(4);
	colors->SetName("EdgeColors");
	// Setup two colors 
	unsigned char red[4] = {255, 0, 0,255};
	unsigned char yellow[4] = {255, 255, 0,0};
	
	// Setup the collision tag array
	vtkSmartPointer<vtkUnsignedCharArray> collisionTag = vtkSmartPointer<vtkUnsignedCharArray>::New();
	collisionTag->SetNumberOfComponents(1);
	collisionTag->SetName("CollisionTag");

	//for (OmMesh::EdgeIter e_it=mesh.edges_begin(); e_it!=mesh.edges_end(); ++e_it) 
	int num_edges = mesh.n_edges();
#pragma omp parallel for	
	for (int i = 0 ; i < num_edges; ++i)
	{
		//OmMesh::EdgeHandle eh = *e_it;
		OmMesh::EdgeHandle eh = mesh.edge_handle(i);	
		
		OmMesh::HalfedgeHandle heh = mesh.halfedge_handle(eh,0);
		if (!heh.is_valid())
			heh = mesh.halfedge_handle(eh,1);
		OmMesh::HalfedgeHandle heh_prev = mesh.prev_halfedge_handle(heh);
		int a = mesh.to_vertex_handle(heh).idx();
		int b = mesh.to_vertex_handle(heh_prev).idx();
		#pragma omp critical
		{
			_graph->AddEdge(a,b);
			colors->InsertNextTupleValue(yellow);
			collisionTag->InsertNextValue(0);
		}
	}
	_graph->GetEdgeData()->AddArray(colors);
	_graph->GetEdgeData()->AddArray(collisionTag);
	this->graph = _graph;
	
	//TagSurroundGraphVertexEdges(_collision_graph,1);
	TagSurroundGraphFaceEdges(_collision_graph);
	clock_t stop = clock();
	timeConsumed += (static_cast<double>(stop)-static_cast<double>(start))/CLOCKS_PER_SEC;
	
	return _graph;
}

vtkSmartPointer<vtkMutableUndirectedGraph> GIMmodTruncate::InitOriginal( vtkSmartPointer<vtkPolyData> _polydata, 															
																		 geodesic::GeodesicAlgorithmExact *_geo , int geoSourceVertexID)
{
	timeConsumed = 0;
	clock_t start = clock();
	seedRemoved = false;
	allFaceRemoved =false;
	firstTruncateDone = false;
	shortenRingsDone = false;
	removePathsBetweenBoundariesFromGraph = false;
	mesh.clear();
	vtkPolydata2OpenMesh(_polydata,&mesh );	
	this->geodesicExact = _geo;	
	this->polydata = _polydata;	
	this->original_polydata = vtkSmartPointer<vtkPolyData>::New();
	this->original_polydata->DeepCopy(_polydata);
	this->original_polydata->BuildLinks();
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
	colors->SetNumberOfComponents(4);
	colors->SetName("EdgeColors");
	// Setup two colors 
	unsigned char red[4] = {255, 0, 0,255};
	unsigned char yellow[4] = {255, 255, 0,255};
	unsigned char purple[4] = {255, 0, 255,255};
	
	// Setup the collision tag array
	vtkSmartPointer<vtkUnsignedCharArray> collisionTag = vtkSmartPointer<vtkUnsignedCharArray>::New();
	collisionTag->SetNumberOfComponents(1);
	collisionTag->SetName("CollisionTag");
	
	//for (OmMesh::EdgeIter e_it=mesh.edges_begin(); e_it!=mesh.edges_end(); ++e_it) 
	int num_edges = mesh.n_edges();
#pragma omp parallel for
	for (int i = 0 ; i < num_edges; ++i)
	{
		//OmMesh::EdgeHandle eh = *e_it;
		OmMesh::EdgeHandle eh = mesh.edge_handle(i);
		
		OmMesh::HalfedgeHandle heh = mesh.halfedge_handle(eh,0);
		if (!heh.is_valid())
			heh = mesh.halfedge_handle(eh,1);
		OmMesh::HalfedgeHandle heh_prev = mesh.prev_halfedge_handle(heh);
		int a = mesh.to_vertex_handle(heh).idx();
		int b = mesh.to_vertex_handle(heh_prev).idx();
		#pragma omp critical
		{
			_graph->AddEdge(a,b);
			colors->InsertNextTupleValue(yellow);
			collisionTag->InsertNextValue(0);
		}

	}
	_graph->GetEdgeData()->AddArray(colors);
	_graph->GetEdgeData()->AddArray(collisionTag);
	this->graph = _graph;	
	clock_t stop = clock();
	timeConsumed += (static_cast<double>(stop)-static_cast<double>(start))/CLOCKS_PER_SEC;
	return _graph;
}

void GIMmodTruncate::Process()
{
	if (!seedRemoved)
	{		
		clock_t start = clock();
		RemoveSeed(false);	
		seedRemoved = true;
	
		RemoveEdgesThatAdjacentOnlyOneFace(false);
		allFaceRemoved = true;
	
		TruncateGraph(false);
		firstTruncateDone = true;

		RemoveOriginalBoundariesFromGraph();
		removePathsBetweenBoundariesFromGraph = true;

		ShorthenRings(false);
		shortenRingsDone = true;

		MergeWithOriginalBoundaries(false);
		
		clock_t stop = clock();
		timeConsumed += (static_cast<double>(stop)-static_cast<double>(start))/CLOCKS_PER_SEC;
		graph->Modified();
		mesh.garbage_collection();
		polydata->RemoveDeletedCells();
		polydata->BuildLinks();
	}
	//GetDiskTopologyPolydata();
	
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
	else if (!removePathsBetweenBoundariesFromGraph)
	{
		RemoveOriginalBoundariesFromGraph();
		removePathsBetweenBoundariesFromGraph = true;
	}
	else if (!shortenRingsDone)
	{
		ShorthenRings(true);
		shortenRingsDone = true;
	}
	else
	{
		MergeWithOriginalBoundaries(false);
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
			unsigned char priority = collisionTag->GetValue(edgeid);
			switch (priority)
			{
			case 0:
				candidate_TagP0Edges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
				 break;
			case 1:
				candidate_TagP1Edges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
				break;
			case 2:
				candidate_TagP2Edges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
				break;
			default:
				throw;
			}
			/*
			if (collisionTag->GetValue(edgeid) == 0)
				candidate_nonTagEdges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
			else
				candidate_TagEdges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
				*/
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
		//OmMesh::Point faceCenter;
		//mesh.calc_face_centroid(removeFaceHandle,faceCenter);
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
				unsigned char priority = collisionTag->GetValue(edgeid);
				switch (priority)
				{
				case 0:
					candidate_TagP0Edges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
					 break;
				case 1:
					candidate_TagP1Edges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
					break;
				case 2:
					candidate_TagP2Edges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
					break;
				default:
					throw;
				}
				/*
				if (collisionTag->GetValue(edgeid) == 0)
					candidate_nonTagEdges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
				else
					candidate_TagEdges.insert(std::pair<double,hedge_data>(best_distance,hedata ));					
					*/
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
	/*
	while (candidate_TagEdges.size() > 0 || candidate_nonTagEdges.size() > 0)
	{
		if (candidate_nonTagEdges.size() > 0)
		{		
			hedge_data checkingEdge = candidate_nonTagEdges.begin()->second;
			candidate_nonTagEdges.erase(candidate_nonTagEdges.begin());
			GenerateCandidateEdgesAndDeleteEdgeFace(checkingEdge,collisionTag);
		
		}
		else if (candidate_TagEdges.size() > 0)
		{
			hedge_data checkingEdge = candidate_TagEdges.begin()->second;
			candidate_TagEdges.erase(candidate_TagEdges.begin());
			GenerateCandidateEdgesAndDeleteEdgeFace(checkingEdge,collisionTag);
		}			
	} //end while
	*/

	
	while (candidate_TagP0Edges.size() > 0 || candidate_TagP1Edges.size() > 0 || candidate_TagP2Edges.size() > 0)
	{
		if (candidate_TagP0Edges.size() > 0)
		{		
			hedge_data checkingEdge = candidate_TagP0Edges.begin()->second;
			candidate_TagP0Edges.erase(candidate_TagP0Edges.begin());
			GenerateCandidateEdgesAndDeleteEdgeFace(checkingEdge,collisionTag);
		
		}
		else if (candidate_TagP1Edges.size() > 0)
		{
			hedge_data checkingEdge = candidate_TagP1Edges.begin()->second;
			candidate_TagP1Edges.erase(candidate_TagP1Edges.begin());
			GenerateCandidateEdgesAndDeleteEdgeFace(checkingEdge,collisionTag);
		}
		else if (candidate_TagP2Edges.size() > 0)
		{
			hedge_data checkingEdge = candidate_TagP2Edges.begin()->second;
			candidate_TagP2Edges.erase(candidate_TagP2Edges.begin());
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

void GIMmodTruncate::TagSurroundGraphFaceEdges(vtkSmartPointer<vtkMutableUndirectedGraph> collision_graph)
{
	
	unsigned char orange[3] = {255, 127, 0};
	unsigned char red[4] = {255, 0, 0,255};
	unsigned char yellow[4] = {255, 255, 0,255};
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
		vtkIdType cedgeid =  graph->GetEdgeId(vid[0],vid[1]);
		collisionTag->SetValue(cedgeid,2);
		edgecolors->SetTupleValue(cedgeid,red);
		std::vector<vtkIdType> adjacent_cellIDs;
		GetNeighborCell(polydata,vid[0],1,adjacent_cellIDs);
		int faceDetected = 0;
		for (int i = 0; i < adjacent_cellIDs.size(); i++)
		{
			if (!cellCheck[adjacent_cellIDs[i]])
			{
				vtkIdType npts;
				vtkIdType *pts;					
				polydata->GetCellPoints(adjacent_cellIDs[i],npts,pts);
				if (npts != 3)
					throw;
				int score = 0;
				for (int k = 0 ; k < 3 ;k++)
				{
					if (pts[k] == vid[0])
						score++;
					if (pts[k] == vid[1])
						score++;
				}

				if (score == 2)
				{
					for (int k = 0 ; k < 3 ;k++)
					{
						vtkIdType edgeid =  graph->GetEdgeId(pts[k%3],pts[(k+1)%3]);
						if (edgeid < 0) 
							throw;
						if (collisionTag->GetValue(edgeid) == 0)
						{
							collisionTag->SetValue(edgeid,1);
							edgecolors->SetTupleValue(edgeid,yellow);
						}
					}
					cellCheck[adjacent_cellIDs[i]] = true;
				}

			}
		}		
	}
}


void GIMmodTruncate::TagSurroundGraphVertexEdges(vtkSmartPointer<vtkMutableUndirectedGraph> collision_graph,int level)
{
	if (level < 1)
		throw;
	unsigned char red[3] = {255, 0, 0};
	unsigned char orange[3] = {255, 127, 0};
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
		vtkIdType cedgeid =  graph->GetEdgeId(vid[0],vid[1]);
		collisionTag->SetValue(cedgeid,2);
		edgecolors->SetTupleValue(cedgeid,red);
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
						if (collisionTag->GetValue(edgeid) == 0)
						{
							collisionTag->SetValue(edgeid,1);
							edgecolors->SetTupleValue(edgeid,orange);
						}
						
					}

					cellCheck[add_cellIDs[j]] = true;
				}
			}
		}		
	}

}

void GIMmodTruncate::ShorthenRings(bool step)
{
	vtkSmartPointer<vtkEdgeListIterator> eit =vtkSmartPointer<vtkEdgeListIterator>::New();
	graph->GetEdges(eit);
	bool *pointCheck = new bool[graph->GetNumberOfVertices()];
	memset(pointCheck,0,sizeof(bool)*graph->GetNumberOfVertices());
	std::vector<edges_path> edgepaths;
	while (eit->HasNext())
	{
		vtkEdgeType e = eit->Next();
		vtkIdType vid[2] = {e.Source,e.Target};
		if (pointCheck[vid[0]] && pointCheck[vid[1]])
			continue;

		if (graph->GetDegree(vid[0]) == 2 && graph->GetDegree(vid[1]) == 2 ) 
		{
			if (boundaryPoints->IsId(vid[0]) >= 0 && boundaryPoints->IsId(vid[1]) >= 0)
			{
				
				continue;
			}
			//find start and end vertices of this edges 
			vtkIdType endofpathvid[2] ;
			vtkIdType prev_endofpathvid[2];
			for (int i = 0 ; i < 2 ;i++)
			{
				vtkIdType nextvid = vid[i];
				vtkIdType prevvid = vid[(i+1)%2];
				
				while (graph->GetDegree(nextvid) == 2)
				{
					pointCheck[nextvid] = true;
					pointCheck[prevvid] = true;
					vtkSmartPointer<vtkAdjacentVertexIterator> subadjacent =vtkSmartPointer<vtkAdjacentVertexIterator>::New();
					graph->GetAdjacentVertices(nextvid,subadjacent);
					vtkIdType nextnextvid = subadjacent->Next();
					if (nextnextvid == prevvid)
						nextnextvid = subadjacent->Next();

					prevvid = nextvid;
					nextvid = nextnextvid;
				}			
				endofpathvid[i] = nextvid;
				prev_endofpathvid[i] = prevvid;
				
			}
			edgepaths.push_back(edges_path(endofpathvid,prev_endofpathvid));
		}
	}

	

	for (int i = 0; i < edgepaths.size(); i++)
	{
		edges_path &path = edgepaths[i];
		vtkIdType remEdgeID =  graph->GetEdgeId(path.endPointID[0],path.pairEdgeEndPointID[0]);
		if (remEdgeID < 0)
			throw;
		else
			graph->RemoveEdge(remEdgeID);

		TruncateGraph(step);
		//return;
		vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra = vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
		dijkstra->SetInputData(original_polydata);
		dijkstra->RepelPathFromVerticesOn();
		
		vtkSmartPointer<vtkPoints> repelPoints = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkEdgeListIterator> eit =vtkSmartPointer<vtkEdgeListIterator>::New();
		graph->GetEdges(eit);
		memset(pointCheck,0,sizeof(bool)*graph->GetNumberOfVertices());
		while (eit->HasNext())
		{
			vtkEdgeType e = eit->Next();
			vtkIdType vid[2] = {e.Source,e.Target};
			for (int j = 0 ; j < 2;j++)
			if (!pointCheck[vid[j]] )
			{
				pointCheck[vid[j]] = true;
				repelPoints->InsertNextPoint( original_polydata->GetPoint(vid[j]));
			}
		}
		dijkstra->SetRepelVertices(repelPoints);
		dijkstra->SetStartVertex(path.pairEdgeEndPointID[0]);
		dijkstra->SetEndVertex(path.pairEdgeEndPointID[1]);
		
		dijkstra->Update();
		

		

		vtkIdList*idlist = dijkstra->GetIdList();
		int numid = idlist->GetNumberOfIds();		
		for (int j = 0 ; j < numid - 1; j++)
		{			
			graph->AddEdge( idlist->GetId(j),idlist->GetId(j+1));
		}
		graph->AddEdge( path.endPointID[0],path.pairEdgeEndPointID[0]);
		graph->AddEdge( path.endPointID[1],path.pairEdgeEndPointID[1]);		
	}

	for (int v = 0; v < graph->GetNumberOfVertices(); v++)
	{
		if (graph->GetDegree(v) == 2 && boundaryPoints->IsId(v) ==-1 )
		{				
			vtkSmartPointer<vtkAdjacentVertexIterator> adjacent =vtkSmartPointer<vtkAdjacentVertexIterator>::New();
			graph->GetAdjacentVertices(v,adjacent);

			vtkIdType ptid0 = adjacent->Next();
			vtkIdType ptid1 = adjacent->Next();
			vtkIdType cellid = GetCellID(original_polydata, v,ptid0,ptid1);
			if (cellid >=0)
			{
				vtkIdType edge0 = graph->GetEdgeId(v,ptid0);
				graph->RemoveEdge(edge0);
					
				vtkIdType edge1 = graph->GetEdgeId(v,ptid1);
				graph->RemoveEdge(edge1);

				graph->AddEdge(ptid0,ptid1);
					
			}
		}
	}
	

	//find nearest vertex in graph to pair edge endpoint
	/*
	for (int j = 0 ; j < 2; j++)
	{
		vtkSmartPointer<vtkIdList> idlist = GetConnectedVertices(original_polydata,path.pairEdgeEndPointID[j]);
			
		vtkIdType bestID;
		double bestDistance = DBL_MAX;
		for (int k = 0; k < idlist->GetNumberOfIds(); k++)
		{
			vtkIdType ptid = idlist->GetId(k);
			if (graph->GetDegree(ptid) > 0)
			{
				double distance = vtkMath::Distance2BetweenPoints(	original_polydata->GetPoint(path.pairEdgeEndPointID[j]), 
																	original_polydata->GetPoint(ptid));
				if (bestDistance > distance)
				{
					bestID = ptid;
					bestDistance = distance;

				}
			}
		}
		//graph->AddEdge( bestID,path.pairEdgeEndPointID[j]);

	}
	*/
	
}


void GIMmodTruncate::MergeWithOriginalBoundaries(bool step)
{
	// Convert the graph to a polydata
	vtkSmartPointer<vtkGraphToPolyData> graphToPolyData = vtkSmartPointer<vtkGraphToPolyData>::New();
	graphToPolyData->SetInputData(graph);
	graphToPolyData->Update();
	
	
	vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter =  vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	connectivityFilter->SetInputData(graphToPolyData->GetOutput());
	connectivityFilter->SetExtractionModeToAllRegions();
	connectivityFilter->Update();

	int numRegions =  connectivityFilter->GetNumberOfExtractedRegions();


	if (numRegions>1)
	{
		//find biggest connectivity.... as main cut path
		connectivityFilter->SetExtractionModeToLargestRegion();
		connectivityFilter->Update();
		vtkSmartPointer<vtkPolyData> mainCutGraph = connectivityFilter->GetOutput();
		
		vtkSmartPointer<vtkDijkstraGraphGeodesicPathMultiStartEndPoints> shortestPath = vtkSmartPointer<vtkDijkstraGraphGeodesicPathMultiStartEndPoints>::New();

	}
	return;
	/*

	vtkSmartPointer<vtkPolyData> process_borderPolydata = vtkSmartPointer<vtkPolyData>::New();
	process_borderPolydata->DeepCopy(borderPolydata);

	vtkLine* _edge = vtkLine::SafeDownCast(process_borderPolydata->GetCell(0));

	
	std::set<int> BorderPoints;	
	std::vector<int> 
	vtkLine* _edge = vtkLine::SafeDownCast(borderEdges->GetOutput()->GetCell(0));
	
	
		for(vtkIdType i = 1; i < numBorderEdges; i++)
		{			
			vtkLine* line = vtkLine::SafeDownCast(borderEdges->GetOutput()->GetCell(i));		
			int numID = line->GetPointIds()->GetNumberOfIds();
			vtkIdType id0 = line->GetPointIds()->GetId(0);
			vtkIdType id1 = line->GetPointIds()->GetId(1);
			vtkIdType oid0 = borderEdges->GetOldIdFromCurrentID(id0);
			vtkIdType oid1 = borderEdges->GetOldIdFromCurrentID(id1);
			BorderPoints.insert(oid0);
			BorderPoints.insert(oid1);		
		}	
		
	*/
	//analyze loops in bouundary edges

	
	//connect each loop to either graph or another loop, if its shortest path is connecting to graph then merge to graph.

}




vtkSmartPointer<vtkPolyData> GIMmodTruncate::GetDiskTopologyPolydata()
{	
	if (!shortenRingsDone)
		return NULL;
	mesh.clear();
	vtkPolydata2OpenMesh(original_polydata,&mesh );	

	std::vector<OmMesh::HalfedgeHandle> cutpath;
	std::map<vtkIdType, int> newpointGen;
	vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
	output->DeepCopy(original_polydata);
	output->BuildLinks();
	
	//create cut path list
	vtkSmartPointer<vtkEdgeListIterator> eit =vtkSmartPointer<vtkEdgeListIterator>::New();
	graph->GetEdges(eit);
	vtkEdgeType e =  eit->Next();
	OmMesh::HalfedgeHandle first_heh =  mesh.find_halfedge(mesh.vertex_handle(e.Source),mesh.vertex_handle(e.Target));
	OmMesh::HalfedgeHandle heh = first_heh;	
	int fromId = mesh.from_vertex_handle(heh).idx();
	int toId = mesh.to_vertex_handle(heh).idx();
	do
	{
		cutpath.push_back(heh);
		newpointGen.insert(std::pair<vtkIdType, int>(toId, graph->GetDegree(toId) - 1));
		//search next halfedge in graph.
		OmMesh::HalfedgeHandle start_heh = mesh.next_halfedge_handle(heh);
		OmMesh::HalfedgeHandle next_heh = start_heh;
		bool found = false;
		do
		{			
			fromId = mesh.from_vertex_handle(next_heh).idx();
			toId = mesh.to_vertex_handle(next_heh).idx();
			if (graph->GetEdgeId(fromId,toId) >= 0)
			{
				heh = next_heh;
				found = true;
				break;
			}
			next_heh = mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(next_heh));
		}
		while(next_heh != start_heh);

		if (!found)
			throw;
	}
	while (heh != first_heh);

	int numEdge = graph->GetNumberOfEdges();
	int numHalfEdge = (int )cutpath.size();
	for (int i = 0; i < numHalfEdge; i++)
	{
		OmMesh::HalfedgeHandle thisHE = cutpath[i];
		OmMesh::HalfedgeHandle nextHE = cutpath[(i+1)%numHalfEdge];
		OmMesh::FaceHandle     nextFH = mesh.face_handle(nextHE);
		fromId = mesh.from_vertex_handle(thisHE).idx();
		toId = mesh.to_vertex_handle(thisHE).idx();
	
		vtkSmartPointer<vtkIdList> cellIDlist = vtkSmartPointer<vtkIdList>::New();
		output->GetCellEdgeNeighbors(	nextFH.idx(), 
										mesh.from_vertex_handle(nextHE).idx(),
										mesh.to_vertex_handle(nextHE).idx(),
										cellIDlist);

		int &degree = newpointGen.at(toId);
		/*
		if ( cellIDlist->GetNumberOfIds() == 0)
			continue;
		*/
		if (degree == 0)
			continue;
		else
			degree--;
		OmMesh::VertexHandle oriVH = mesh.to_vertex_handle(thisHE);
		OmMesh::VertexHandle newVH = mesh.add_vertex( mesh.point(oriVH));
		vtkIdType oldVertexID = oriVH.idx();
		vtkIdType newVertexID = output->GetPoints()->InsertNextPoint( output->GetPoint(oldVertexID));

		
		OmMesh::HalfedgeHandle start_heh = mesh.next_halfedge_handle(thisHE);
		OmMesh::HalfedgeHandle next_heh = start_heh;
		bool found = false;
		do
		{			
			int faceId = mesh.face_handle(next_heh).idx();				
			output->ReplaceCellPoint( faceId,oldVertexID, newVertexID);		
			
			if (next_heh == nextHE)
			{
				found = true;
				break;
			}
			next_heh = mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(next_heh));
		}
		while(next_heh != start_heh);

		if (!found)
			throw;
		output->BuildLinks();


		
	}
	
	/*
	for (int i = 0 ; i < num_vertex_cut_route; i++)
	{
		
		IDCutHedge::AppendVF(vertex_route[i],face_route[i],0,pCutHedgeT);
		//Search for twin
		
		IDCutHedge *justCreated = (IDCutHedge *)pCutHedgeT->back;		
		//Search for duplicate ID		
		search = justCreated;
		while(back(search) != pCutHedgeH)
		{
			search = back(search);

			if (search->ID == justCreated->ID)
			{
				justCreated->DuplicateDegree = ((IDCutHedge *)search)->DuplicateDegree + 1;
				m_num_boundarySurfacePolarVertex++;;//num_boundarySurfacePolarVertex[degree_count]++;
				numDuplicateVertices++;
				break;
			}
		}

		if ( back(pCutHedgeT)->FaceID == back(back(pCutHedgeT))->FaceID)
		{
			m_numValen2BoundaryPoint++;
		}
	}

	*/
	return output;
}


void GIMmodTruncate::RemoveOriginalBoundariesFromGraph()
{
	//find border points/edges
	vtkSmartPointer<vtkFeatureEdgesEx> borderEdges = vtkSmartPointer<vtkFeatureEdgesEx>::New();
	borderEdges->SetInputData(original_polydata);
	borderEdges->FeatureEdgesOff();
	borderEdges->BoundaryEdgesOn();
	borderEdges->ManifoldEdgesOff();
	borderEdges->NonManifoldEdgesOff();
	borderEdges->Update();	
	
	vtkSmartPointer<vtkPolyData> borderPolydata = borderEdges->GetOutput();
	vtkSmartPointer<vtkCellArray> lines= borderPolydata->GetLines();	
	vtkSmartPointer<vtkPoints> points = borderPolydata->GetPoints();
	vtkIdList *borderPointID = borderEdges->GetOldIdList();
	boundaryPoints = vtkSmartPointer<vtkIdList>::New();
	boundaryPoints->DeepCopy(borderPointID);
	int numBorderEdges = lines->GetNumberOfCells();
	int numBorderVertices = points->GetNumberOfPoints();
	if (numBorderEdges == 0) //if no original boundary .... exit		
		return;
	for(vtkIdType i = 0; i < numBorderVertices; i++)
	{			
		vtkIdType old_id = borderEdges->GetOldIdFromCurrentID(i);
		int inDegree =  graph->GetInDegree(old_id) ;
		if (inDegree >= 3)
		{
			//search for end point is not boundary points
			for (int iter_i = 0; iter_i < inDegree ; iter_i++)
			{
				vtkInEdgeType in_e = graph->GetInEdge(old_id, iter_i);
				if (borderPointID->IsId(in_e.Source) == -1)
				{
					//end point is not boundary point

					//temp remove that edge
					int numGraphPoint = graph->GetNumberOfVertices();
					graph->RemoveEdge(in_e.Id);
					

					vtkSmartPointer<vtkBoostBreadthFirstSearch> bfs = vtkSmartPointer<vtkBoostBreadthFirstSearch>::New();
					bfs->SetInputData(graph);				
					bfs->SetOriginVertex(old_id);
					bfs->Update();
					vtkIntArray *bfsArray = vtkIntArray::SafeDownCast(bfs->GetOutput()->GetVertexData()->GetArray("BFS"));
					
					int value = bfsArray->GetValue(in_e.Source);
					if (value == VTK_INT_MAX)
					{
						/*
						//unconnected , safe to remove
						//remove a boundary edge around this point
						int _inDegree =  graph->GetInDegree(old_id);
						for (int _iter_i = 0; _iter_i < _inDegree ; _iter_i++)
						{
							vtkInEdgeType _in_e = graph->GetInEdge(old_id, _iter_i);
							if (borderPointID->IsId(_in_e.Source) >= 0)
							{
								graph->RemoveEdge(_in_e.Id);
								break;
							}
						}
						*/
						break;
					}
					else
					{
						//not safe to remove 
						//restore edge
						graph->AddEdge(old_id,in_e.Source);
					}
					
					


				}
				else if (graph->GetInDegree(in_e.Source) >= 3)
				{
					//if end point is boundary point  that have valence 3
					//test for remove


					//temp remove that edge
					int numGraphPoint = graph->GetNumberOfVertices();
					graph->RemoveEdge(in_e.Id);
					

					vtkSmartPointer<vtkBoostBreadthFirstSearch> bfs = vtkSmartPointer<vtkBoostBreadthFirstSearch>::New();
					bfs->SetInputData(graph);				
					bfs->SetOriginVertex(old_id);
					bfs->Update();
					vtkIntArray *bfsArray = vtkIntArray::SafeDownCast(bfs->GetOutput()->GetVertexData()->GetArray("BFS"));
					
					int value = bfsArray->GetValue(in_e.Source);
					if (value == VTK_INT_MAX)
					{
						//unconnected , safe to remove
						//remove a boundary edge around this point
						/*
						int _inDegree =  graph->GetInDegree(old_id);
						for (int _iter_i = 0; _iter_i < _inDegree ; _iter_i++)
						{
							vtkInEdgeType _in_e = graph->GetInEdge(old_id, _iter_i);
							if (borderPointID->IsId(_in_e.Source) >= 0)
							{
								graph->RemoveEdge(_in_e.Id);
								break;
							}
						}
						*/
						break;
					}
					else
					{
						//not safe to remove 
						//restore edge
						graph->AddEdge(old_id,in_e.Source);
					}

				}

			}
		}
	}

	TruncateGraph(false);
}