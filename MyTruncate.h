#pragma once
#include "VTK_Header.h"
#include "OpenMesh_Header.h"
#include "geodesic\geodesic_algorithm_exact.h"



class MyTruncate
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
	MyTruncate(void);
	virtual ~MyTruncate(void);
	vtkSmartPointer<vtkMutableUndirectedGraph> Init( vtkSmartPointer<vtkPolyData> polydata, geodesic::GeodesicAlgorithmExact *geo);
	void Step();
	void Process();
	vtkSmartPointer<vtkMutableUndirectedGraph> GetGraph();
protected:
	vtkWeakPointer<vtkPolyData> polydata;
	geodesic::GeodesicAlgorithmExact *geodesicExact;
	OmMesh mesh;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph;
	bool seedRemoved;
	bool allFaceRemoved;
	bool firstTruncateDone;
	bool branchRemoved;
	bool largestGraphDone;
	std::map<double,hedge_data> candidate_edges;

	void RemoveSeed();
	void RemoveBoundarySeed(bool step);
	void RemoveEdgesThatAdjacentOnlyOneFace(bool step);
	void EliminateSelfCycle(bool step);
	void EliminateBranchPath(bool step);
	void TruncateGraph(bool step);
	void EliminateUnusedPath(bool step);
};

