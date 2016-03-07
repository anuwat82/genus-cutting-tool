#pragma once
#include "geodesic\geodesic_algorithm_exact.h"
#include "VTK_Header.h"
#include "GIMmodTruncate.h"

//extern MyTruncate gimTruncate;
extern GIMmodTruncate modTruncate;
extern GIMmodTruncate originalTruncate;
extern vtkWeakPointer<vtkPolyData> polydata;
extern vtkWeakPointer<vtkPoints> source_point;
extern vtkSmartPointer<vtkPolyData> disk_polydata;
extern vtkWeakPointer<vtkPolyDataMapper> mapper;

extern vtkWeakPointer<vtkActor> actorMainPoly;
extern vtkWeakPointer<vtkActor> actorEdge1;
extern vtkWeakPointer<vtkActor> actorEdge2;
extern vtkWeakPointer<vtkActor> actorEdge3;
extern vtkWeakPointer<vtkActor> actorEdge4;
extern vtkWeakPointer<vtkActor> actorEdge5;

extern vtkWeakPointer<vtkActor> actorPoly1;
extern vtkWeakPointer<vtkActor> actorPoly2;
extern vtkWeakPointer<vtkActor> actorPoly3;
extern geodesic::Mesh geosesic_mesh;
extern geodesic::GeodesicAlgorithmExact *exact_algorithm;
extern int sourceVertex;


struct collisionEdgeInfo
{
	int v0;
	int v1;
	int type;
};

void ColorMeshVertice(vtkDoubleArray *scalar);
void ColorMeshFace(vtkDoubleArray *scalar);
void ColoredPoint(vtkSmartPointer<vtkRenderer> _renderer , double pts[3],  double r, double g, double b);

void LoadToGeodesic(vtkSmartPointer<vtkPolyData> polydata);
bool GenerateGeodesicDistance(geodesic::GeodesicAlgorithmExact &geo ,int vertexID , std::vector<collisionEdgeInfo>& collision_edges);
vtkSmartPointer<vtkMutableUndirectedGraph> GenerateCollisionEdgeGraph(vtkSmartPointer<vtkPolyData> source_polydata, std::vector<collisionEdgeInfo>& collision_edges);
vtkSmartPointer<vtkMutableUndirectedGraph> TruncateGraph(vtkSmartPointer<vtkMutableUndirectedGraph> source_graph);
vtkSmartPointer<vtkMutableUndirectedGraph> StraightupGraph(vtkSmartPointer<vtkMutableUndirectedGraph> source_graph, OmMesh& mesh);
vtkSmartPointer<vtkMutableUndirectedGraph> CleanGraph(vtkSmartPointer<vtkMutableUndirectedGraph> source_graph, OmMesh& mesh);
vtkSmartPointer<vtkPolyData> CreateSurroundGraphPolydata(vtkSmartPointer<vtkMutableUndirectedGraph> source_graph, vtkSmartPointer<vtkPolyData> source_polydata,int level = 1);
vtkSmartPointer<vtkMutableUndirectedGraph> GIMtruncate(vtkSmartPointer<vtkPolyData> polydata, geodesic::GeodesicAlgorithmExact &geo);
vtkSmartPointer<vtkActor> CreateBeforeTruncatePipeline(vtkSmartPointer<vtkMutableUndirectedGraph> BTGraph,std::vector<collisionEdgeInfo>& collision_edges);
vtkSmartPointer<vtkActor> CreateAfterTruncatePipeline(vtkSmartPointer<vtkMutableUndirectedGraph> ATGraph,std::vector<collisionEdgeInfo>& collision_edges);
vtkSmartPointer<vtkActor> CreateStrightUpPipeline(vtkSmartPointer<vtkMutableUndirectedGraph> SUGraph);
vtkSmartPointer<vtkActor> CreateCleanPipeline(vtkSmartPointer<vtkMutableUndirectedGraph> CleanGraph);
vtkSmartPointer<vtkActor> CreateFinalPipeline(vtkSmartPointer<vtkMutableUndirectedGraph> FinalGraph);
vtkSmartPointer<vtkMutableUndirectedGraph> CreateBoundaryGraph(vtkSmartPointer<vtkPolyData> polydata);

void InitialGeodesic(vtkSmartPointer<vtkPolyData> polydata ,int sourceVertexID);
void ScreenShot(vtkRenderWindow* renderWindow);
void AskForSaveSqp(vtkFloatArray *texCoord);
void AskForSaveParameterizationPLY(vtkPolyData* diskTopology, vtkFloatArray *texCoord);
vtkSmartPointer<vtkPolyData> CleanForUnrefVertex(vtkSmartPointer<vtkPolyData> input);

float GetL2ErrorMS(vtkPolyData* diskTopology, vtkFloatArray *texCoord);
