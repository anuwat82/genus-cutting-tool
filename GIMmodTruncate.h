#pragma once
#include "VTK_Header.h"
#include "OpenMesh_Header.h"
#include "geodesic\geodesic_algorithm_exact.h"



class GIMmodTruncate
{
public:
	struct hedge_data
	{	
		vtkIdType startID;
		vtkIdType endID;
		hedge_data( vtkIdType start, vtkIdType end)
		{		
			startID = start;
			endID = end;
		}
	};

	struct edges_path
	{	
		vtkIdType endPointID[2];
		vtkIdType pairEdgeEndPointID[2]; //
		edges_path(){}
		edges_path( vtkIdType epid[2], vtkIdType peepid[2])
		{
			endPointID[0] = epid[0];
			endPointID[1] = epid[1];
			pairEdgeEndPointID[0] = peepid[0];
			pairEdgeEndPointID[1] = peepid[1];
		}
		
	};
	GIMmodTruncate(void);
	virtual ~GIMmodTruncate(void);
	vtkSmartPointer<vtkMutableUndirectedGraph> Init( vtkSmartPointer<vtkPolyData> polydata, 
													 vtkSmartPointer<vtkMutableUndirectedGraph> collision_graph,
													 geodesic::GeodesicAlgorithmExact *geo , int geoSourceVertexID);
	vtkSmartPointer<vtkMutableUndirectedGraph> InitOriginal( vtkSmartPointer<vtkPolyData> polydata, 															 
															 geodesic::GeodesicAlgorithmExact *geo , int geoSourceVertexID);
	void Step();
	void Process();
	vtkSmartPointer<vtkMutableUndirectedGraph> GetGraph();
	vtkSmartPointer<vtkPolyData> GetDiskTopologyPolydata();

	bool isReadyToCut();
protected:
	vtkWeakPointer<vtkPolyData> polydata;
	vtkSmartPointer<vtkPolyData> original_polydata;
	geodesic::GeodesicAlgorithmExact *geodesicExact;
	OmMesh mesh;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph;
	int geodesicSourceVertexID;
	bool seedRemoved;
	bool allFaceRemoved;
	bool firstTruncateDone;
	bool shortenRingsDone;
	bool largestGraphDone;
	std::multimap<double,hedge_data> candidate_nonTagEdges;
	std::multimap<double,hedge_data> candidate_TagEdges;

	
	std::multimap<double,hedge_data> candidate_TagP0Edges;
	std::multimap<double,hedge_data> candidate_TagP2Edges;
	std::multimap<double,hedge_data> candidate_TagP1Edges;

	void RemoveSeed(bool ste);
	//void RemoveBoundarySeed(bool step);
	void RemoveEdgesThatAdjacentOnlyOneFace(bool step);
	void GenerateCandidateEdgesAndDeleteEdgeFace(hedge_data &checkingEdge ,vtkSmartPointer<vtkUnsignedCharArray> collisionTag);
	void EliminateSelfCycle(bool step);
	void EliminateBranchPath(bool step);
	void TruncateGraph(bool step);
	void EliminateUnusedPath(bool step);

	void ShorthenRings(bool step);

	void TagSurroundGraphFaceEdges(vtkSmartPointer<vtkMutableUndirectedGraph> collision_graph);// tag surrond edges of face that adjacent on an edge in graph
	void TagSurroundGraphVertexEdges(vtkSmartPointer<vtkMutableUndirectedGraph> collision_graph,int level);  // tag surrond edges of face around a vertex in graph
	void GetNeighborCell(vtkPolyData* source_polydata, vtkIdType vid , int level, std::vector<vtkIdType> &storeIDs);
};

