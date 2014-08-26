#pragma once
#include "VTK_Header.h"
#include "OpenMesh_Header.h"
#include "geodesic\geodesic_algorithm_exact.h"



class GIMTruncate
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
	GIMTruncate(void);
	virtual ~GIMTruncate(void);
	vtkSmartPointer<vtkMutableUndirectedGraph> Init( vtkSmartPointer<vtkPolyData> polydata, geodesic::GeodesicAlgorithmExact *geo);
	void Step();
	vtkSmartPointer<vtkMutableUndirectedGraph> GetGraph();
protected:
	vtkWeakPointer<vtkPolyData> polydata;
	geodesic::GeodesicAlgorithmExact *geodesicExact;
	OmMesh mesh;
	vtkSmartPointer<vtkMutableUndirectedGraph> graph;
	bool seedRemoved;
	bool allFaceRemoved;
	bool firstTruncateDone;
	std::map<double,hedge_data> candidate_edges;

	void RemoveSeed();
	void RemoveBoundarySeed();
	void EliminateSelfCycle();
	void EliminateBranchPath();
	void TruncateGraph();
};

