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
	GIMmodTruncate(void);
	virtual ~GIMmodTruncate(void);
	vtkSmartPointer<vtkMutableUndirectedGraph> Init( vtkSmartPointer<vtkPolyData> polydata, 
													 vtkSmartPointer<vtkMutableUndirectedGraph> collision_graph,
													 geodesic::GeodesicAlgorithmExact *geo , int geoSourceVertexID);
	void Step();
	void Process();
	vtkSmartPointer<vtkMutableUndirectedGraph> GetGraph();
protected:
	vtkWeakPointer<vtkPolyData> polydata;
	geodesic::GeodesicAlgorithmExact *geodesicExact;
	OmMesh mesh;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph;
	int geodesicSourceVertexID;
	bool seedRemoved;
	bool allFaceRemoved;
	bool firstTruncateDone;
	bool branchRemoved;
	bool largestGraphDone;
	std::map<double,hedge_data> candidate_nonTagEdges;
	std::map<double,hedge_data> candidate_TagEdges;

	void RemoveSeed(bool ste);
	//void RemoveBoundarySeed(bool step);
	void RemoveEdgesThatAdjacentOnlyOneFace(bool step);
	void GenerateCandidateEdgesAndDeleteEdgeFace(hedge_data &checkingEdge ,vtkSmartPointer<vtkUnsignedCharArray> collisionTag);
	void EliminateSelfCycle(bool step);
	void EliminateBranchPath(bool step);
	void TruncateGraph(bool step);
	void EliminateUnusedPath(bool step);

	void TagSurroundGraphEdges(vtkSmartPointer<vtkMutableUndirectedGraph> collision_graph,int level);
	void GetNeighborCell(vtkPolyData* source_polydata, vtkIdType vid , int level, std::vector<vtkIdType> &storeIDs);
};

