
#include "stdafx.h"

#include <Windows.h>
#include "utils.h"
#include "VTK_Header.h"
#include "MouseInteractorStylePP.h"
#include "vtkOFFReader.h"
#include <ostream>
#include <string>     // std::string, std::stoi
#include "geodesic\geodesic_algorithm_exact.h"
#include "MyTruncate.h"
#include "GIMmodTruncate.h"
#include "Parameterization\PolygonsData.h"
MyTruncate gimTruncate;
GIMmodTruncate modTruncate;
GIMmodTruncate originalTruncate;
vtkWeakPointer<vtkPolyData> polydata;
vtkWeakPointer<vtkPoints> source_point;

vtkWeakPointer<vtkActor> actorMainPoly;
vtkWeakPointer<vtkActor> actorEdge1;
vtkWeakPointer<vtkActor> actorEdge2;
vtkWeakPointer<vtkActor> actorEdge3;
vtkWeakPointer<vtkActor> actorEdge4;
vtkWeakPointer<vtkActor> actorEdge5;

vtkWeakPointer<vtkActor> actorPoly1;
vtkWeakPointer<vtkActor> actorPoly2;
vtkWeakPointer<vtkActor> actorPoly3;
geodesic::Mesh geosesic_mesh;
geodesic::GeodesicAlgorithmExact *exact_algorithm;
struct collisionEdgeInfo
{
	int v0;
	int v1;
	int type;
};
void ColoredPoint(vtkSmartPointer<vtkRenderer> _renderer , double pts[3],  double r, double g, double b);
void keyPressCallbackFunc(vtkObject*, unsigned long eid, void* clientdata, void *calldata);

void pickCallbackFunc(vtkObject*, unsigned long eid, void* clientdata, void *calldata);
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

void Process(vtkSmartPointer<vtkPolyData> polydata ,int sourceVertexID);
void ScreenShot(vtkRenderWindow* renderWindow);
std::string GetFileExtension(const std::string& FileName)
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}
/*
void Process(vtkSmartPointer<vtkPolyData> polydata , int sourceVertexID )
{
	std::vector<collisionEdgeInfo> collision_edges;
	OmMesh mesh;
	vtkPolydata2OpenMesh(polydata,&mesh );
	if (polydata)	
		LoadToGeodesic(polydata);

	int numVertex = polydata->GetNumberOfPoints();
	if (sourceVertexID >= numVertex)
	{
		std::cout << "invalid source vertex id ... reset to id 0." << endl;
		sourceVertexID = 0;
	}
	std::cout << "==========================" << endl;
	std::cout << "source vertex id " << sourceVertexID << endl;

	geodesic::GeodesicAlgorithmExact exact_algorithm(&geosesic_mesh);
	GenerateGeodesicDistance(exact_algorithm,sourceVertexID,collision_edges);
	vtkSmartPointer<vtkMutableUndirectedGraph> collisionEdgesGraphBeforeTruncate = GenerateCollisionEdgeGraph(polydata,collision_edges);
	vtkSmartPointer<vtkMutableUndirectedGraph> collisionEdgesGraphAfterTruncate  = TruncateGraph(collisionEdgesGraphBeforeTruncate);
	vtkSmartPointer<vtkMutableUndirectedGraph> collisionEdgesStraight			 = StraightupGraph(collisionEdgesGraphAfterTruncate,mesh);
	vtkSmartPointer<vtkMutableUndirectedGraph> collisionEdgesClean				 = CleanGraph(collisionEdgesStraight,mesh);
	


	vtkSmartPointer<vtkMutableUndirectedGraph> graph  = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	graph->DeepCopy(collisionEdgesClean);
	int adjnum = graph->GetDegree(140);
	adjnum = graph->GetDegree(152);
	
	int numGraphEdges = 0;
	do 
	{
		numGraphEdges = graph->GetNumberOfEdges();
		graph = TruncateGraph(graph);
		graph = StraightupGraph(graph,mesh);
		graph = CleanGraph(graph,mesh);
	
	}
	while (numGraphEdges > graph->GetNumberOfEdges());
	
	vtkSmartPointer<vtkActor> edge_actor1 = CreateBeforeTruncatePipeline(collisionEdgesGraphBeforeTruncate, collision_edges); 
	vtkSmartPointer<vtkActor> edge_actor2 = CreateAfterTruncatePipeline(collisionEdgesGraphAfterTruncate, collision_edges);
	vtkSmartPointer<vtkActor> edge_actor3 = CreateStrightUpPipeline(collisionEdgesStraight);
	vtkSmartPointer<vtkActor> edge_actor4 = CreateCleanPipeline(collisionEdgesClean);
	vtkSmartPointer<vtkActor> edge_actor5 = CreateCleanPipeline(graph);
	if (edge_actor1->GetReferenceCount() == 1)
		edge_actor1->SetReferenceCount(2);
	if (edge_actor2->GetReferenceCount() == 1)
		edge_actor2->SetReferenceCount(2);
	if (edge_actor3->GetReferenceCount() == 1)
		edge_actor3->SetReferenceCount(2);
	if (edge_actor4->GetReferenceCount() == 1)
		edge_actor4->SetReferenceCount(2);
	if (edge_actor5->GetReferenceCount() == 1)
		edge_actor5->SetReferenceCount(2);
	
	if (actorEdge1)
		actorEdge1->ShallowCopy(edge_actor1);
	else
		actorEdge1 = edge_actor1;
	if (actorEdge2)
		actorEdge2->ShallowCopy(edge_actor2);
	else
		actorEdge2 = edge_actor2;

	if (actorEdge3)
		actorEdge3->ShallowCopy(edge_actor3);
	else
		actorEdge3 = edge_actor3;

	if (actorEdge4)
		actorEdge4->ShallowCopy(edge_actor4);
	else
		actorEdge4 = edge_actor4;

	if (actorEdge5)
		actorEdge5->ShallowCopy(edge_actor5);
	else
		actorEdge5 = edge_actor5;

	vtkSmartPointer<vtkPolyData> surroundPolydata = CreateSurroundGraphPolydata(collisionEdgesGraphAfterTruncate,polydata);
	GIMtruncate(surroundPolydata,exact_algorithm);
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData (surroundPolydata);
 
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetEdgeVisibility(1);
	actor->GetProperty()->SetLineWidth(0.5);
	if (actor->GetReferenceCount() == 1)
		actor->SetReferenceCount(2);

	if (actorPoly2)
		actorPoly2->ShallowCopy(actor);
	else
		actorPoly2 = actor;
	 

	
}
*/

void Process(vtkSmartPointer<vtkPolyData> polydata , int sourceVertexID )
{
	std::vector<collisionEdgeInfo> collision_edges;
	OmMesh mesh;
	vtkPolydata2OpenMesh(polydata,&mesh );
	if (polydata)	
		LoadToGeodesic(polydata);

	int numVertex = polydata->GetNumberOfPoints();
	if (sourceVertexID >= numVertex)
	{
		std::cout << "invalid source vertex id ... reset to id 0." << endl;
		sourceVertexID = 0;
	}
	std::cout << "==========================" << endl;
	std::cout << "source vertex id " << sourceVertexID << endl;
	
	exact_algorithm = new geodesic::GeodesicAlgorithmExact(&geosesic_mesh);
	GenerateGeodesicDistance(*exact_algorithm,sourceVertexID,collision_edges);
	std::cout << "number of edges:" << exact_algorithm->mesh()->edges().size() << endl; 

	vtkSmartPointer<vtkMutableUndirectedGraph> collisionEdgesGraphBeforeTruncate = GenerateCollisionEdgeGraph(polydata,collision_edges);
	vtkSmartPointer<vtkMutableUndirectedGraph> collisionEdgesGraphAfterTruncate  = TruncateGraph(collisionEdgesGraphBeforeTruncate);
	//vtkSmartPointer<vtkPolyData> surroundPolydata = CreateSurroundGraphPolydata(collisionEdgesGraphAfterTruncate,polydata,2);	
	//vtkSmartPointer<vtkMutableUndirectedGraph> graphState3 = gimTruncate.Init(surroundPolydata,exact_algorithm);
	vtkSmartPointer<vtkPolyData> truncateModPolydata = vtkSmartPointer<vtkPolyData>::New();
	truncateModPolydata->DeepCopy(polydata);
	truncateModPolydata->BuildLinks();
	vtkSmartPointer<vtkPolyData> truncateOriginalPolydata = vtkSmartPointer<vtkPolyData>::New();
	truncateOriginalPolydata->DeepCopy(polydata);
	truncateOriginalPolydata->BuildLinks();
	vtkSmartPointer<vtkMutableUndirectedGraph> graphState3 = modTruncate.Init(truncateModPolydata,collisionEdgesGraphAfterTruncate,exact_algorithm,sourceVertexID);
	vtkSmartPointer<vtkMutableUndirectedGraph> graphState4 = originalTruncate.InitOriginal(truncateOriginalPolydata,exact_algorithm,sourceVertexID);
	
	{
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData (truncateModPolydata);
 
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetEdgeVisibility(1);
		actor->GetProperty()->SetLineWidth(0.5);
		actor->GetProperty()->SetColor(237.0/255.0 , 249.0/255.0,200.0/255.0);
		if (actor->GetReferenceCount() == 1)
			actor->SetReferenceCount(2);
		if (actorPoly3)
			actorPoly3->ShallowCopy(actor);
		else
			actorPoly3 = actor;
	}

	{
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputData (truncateOriginalPolydata);
 
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(mapper);
		actor->GetProperty()->SetEdgeVisibility(1);
		actor->GetProperty()->SetLineWidth(0.5);
		actor->GetProperty()->SetColor(174.0/255.0 , 225.0/255.0,238.0/255.0);
		if (actor->GetReferenceCount() == 1)
			actor->SetReferenceCount(2);
		if (actorPoly2)
			actorPoly2->ShallowCopy(actor);
		else
			actorPoly2 = actor;
	}

	vtkSmartPointer<vtkActor> edge_actor1 = CreateBeforeTruncatePipeline(collisionEdgesGraphBeforeTruncate, collision_edges); 
	vtkSmartPointer<vtkActor> edge_actor2 = CreateAfterTruncatePipeline(collisionEdgesGraphAfterTruncate, collision_edges);
	vtkSmartPointer<vtkActor> edge_actor3 = CreateStrightUpPipeline(graphState3);
	vtkSmartPointer<vtkActor> edge_actor4 = CreateStrightUpPipeline(graphState4);
	
	if (edge_actor1->GetReferenceCount() == 1)
		edge_actor1->SetReferenceCount(2);
	if (edge_actor2->GetReferenceCount() == 1)
		edge_actor2->SetReferenceCount(2);
	if (edge_actor3->GetReferenceCount() == 1)
		edge_actor3->SetReferenceCount(2);	
	if (edge_actor4->GetReferenceCount() == 1)
		edge_actor4->SetReferenceCount(2);
	/*
	if (edge_actor5->GetReferenceCount() == 1)
		edge_actor5->SetReferenceCount(2);
	*/
	if (actorEdge1)
		actorEdge1->ShallowCopy(edge_actor1);
	else
		actorEdge1 = edge_actor1;
	if (actorEdge2)
		actorEdge2->ShallowCopy(edge_actor2);
	else
		actorEdge2 = edge_actor2;

	if (actorEdge3)
		actorEdge3->ShallowCopy(edge_actor3);
	else
		actorEdge3 = edge_actor3;
	
	if (actorEdge4)
		actorEdge4->ShallowCopy(edge_actor4);
	else
		actorEdge4 = edge_actor4;
	/*
	if (actorEdge5)
		actorEdge5->ShallowCopy(edge_actor5);
	else
		actorEdge5 = edge_actor5;

	*/
}


int main(int argc, char* argv[])
{
	int sourceVertex  = 0;
	std::string filename;
	if (argc >= 2)
	{
		filename = std::string(argv[1]);
		if (argc >= 3)
		{		
			sourceVertex = _tstoi(argv[2]);
			if (errno != 0)
			{
				std::cout << "Cannot convert parameter #2 to source vertex id" << endl;
				sourceVertex = 0;
			}
		}
	}
	else
	{
		//open dialog for 
		if (!GetFileName(filename,"PLY Files\0*.ply\0OFF files\0*.off\0All Files\0*.*\0"))
		{
			std::cout << "Exit due to no model file input." << endl;
			return -1;
		}
	}
	std::string ext = GetFileExtension(filename);
	vtkSmartPointer<vtkPolyDataAlgorithm> modelReader;
	if (ext == "ply")
	{
		vtkSmartPointer<vtkPLYReader> PLYReader = vtkSmartPointer<vtkPLYReader>::New();
		PLYReader->SetFileName(filename.c_str());
		PLYReader->Update();
		modelReader = PLYReader;		
	}
	else if (ext == "off")
	{
		vtkSmartPointer<vtkOFFReader> OFFReader = vtkSmartPointer<vtkOFFReader>::New();
		OFFReader->SetFileName(filename.c_str());
		OFFReader->Update();
		modelReader = OFFReader;		
	}
	polydata = modelReader->GetOutput();
	Process(modelReader->GetOutput() , sourceVertex);


	// Create the tree
	vtkSmartPointer<vtkOBBTree> obbTree = vtkSmartPointer<vtkOBBTree>::New();  
	obbTree->SetDataSet(polydata);
	obbTree->BuildLocator();
	
	

	
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(modelReader->GetOutputPort());
 
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actorMainPoly = actor;
	actor->GetProperty()->SetEdgeVisibility(0);
	actor->GetProperty()->SetLineWidth(0.1);
	actor->GetProperty()->SetOpacity(0.75); //0.65
	// Visualize
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->StereoCapableWindowOff();
	renderWindow->StereoRenderOff(); 	
	
	//depth peeling
	renderWindow->SetAlphaBitPlanes(1);
	renderWindow->SetMultiSamples(0);
	renderer->SetUseDepthPeeling(1);
	renderer->SetMaximumNumberOfPeels(100);
	renderer->SetOcclusionRatio(0.1);
	
	
	
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<MouseInteractorStylePP> TrackballStyle = vtkSmartPointer<MouseInteractorStylePP>::New();		
	vtkSmartPointer<vtkPointPicker> picker =  vtkSmartPointer<vtkPointPicker>::New();
	double tol = picker->GetTolerance();
	picker->SetTolerance(tol*0.25);

		
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderWindowInteractor->SetInteractorStyle(TrackballStyle);
	renderWindowInteractor->SetPicker(picker);
	
	//renderWindowInteractor->SetPicker(cellPicker);
	TrackballStyle->SetPickColor(1.0,0.0,0.0);
	actorPoly1 = actor;
	renderer->AddActor(actorPoly1);
	renderer->AddActor(actorEdge2);
	
	/*
	// Initialize the representation
	for (int level = 0; level <= 2; level++)
	{
		vtkSmartPointer<vtkPolyData> obbpolydata = vtkSmartPointer<vtkPolyData>::New();
		obbTree->GenerateRepresentation(level, obbpolydata);
 
		vtkSmartPointer<vtkPolyDataMapper> obbtreeMapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
		obbtreeMapper->SetInputData(obbpolydata); 
	
		vtkSmartPointer<vtkActor> obbtreeActor = vtkSmartPointer<vtkActor>::New();
		obbtreeActor->SetMapper(obbtreeMapper);
		obbtreeActor->GetProperty()->SetInterpolationToFlat();
		obbtreeActor->GetProperty()->SetRepresentationToWireframe();

		double color[3] = {0,0,0};
		color[level] +=0.5;
		obbtreeActor->GetProperty()->SetColor( color);
		renderer->AddActor(obbtreeActor);
	}
	*/
	
	vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	keypressCallback->SetCallback ( keyPressCallbackFunc );
	renderWindowInteractor->AddObserver ( vtkCommand::KeyPressEvent, keypressCallback );
	

	vtkSmartPointer<vtkCallbackCommand> pickCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	pickCallback->SetCallback ( pickCallbackFunc );
	pickCallback->SetClientData(TrackballStyle);
	picker->AddObserver(vtkCommand::EndPickEvent,pickCallback);	
	
	{
		vtkCamera *cam  = renderer->GetActiveCamera();
		double pos[3] = {7.21609, -25.5471, 20.343};
		
		double up[3] = {-0.17565, 0.582402, 0.793697};
		
		cam->SetPosition(pos);
		cam->SetViewUp(up);
		
	}

	ColoredPoint( renderer,modelReader->GetOutput()->GetPoint(sourceVertex), 1.0,0.5,0.0);

	renderWindow->Render(); 	

	cout << "===============================" << endl <<
		"Function Key:" << endl <<
		"F6 = homotopy cutting" << endl <<
		"F7 = iterative augmented cutting" << endl <<
		"F8 = square parameterization (optimal)" << endl << endl <<
		
		"F10 = save screenshot" << endl << 
		"===============================" << endl;


	renderWindowInteractor->Start();
	
	return 0;
}

void keyPressCallbackFunc(vtkObject* caller, unsigned long eid, void* clientdata, void *calldata)
{
	vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
	vtkRenderer *renderer = static_cast<vtkRenderer*>(iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
 
	
	char *ch = iren->GetKeySym();
	switch (*ch)
	{
		case 'r':
			{
				vtkCamera *cam  = renderer->GetActiveCamera();
				cout<< *cam << endl;
			}
			break;
		case '1':
			renderer->RemoveActor(actorEdge1);
			renderer->RemoveActor(actorEdge2);
			renderer->RemoveActor(actorEdge3);
			renderer->RemoveActor(actorEdge4);
			renderer->RemoveActor(actorEdge5);
			renderer->AddActor(actorEdge1);
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			break;
		case '2':
			renderer->RemoveActor(actorEdge1);
			renderer->RemoveActor(actorEdge2);
			renderer->RemoveActor(actorEdge3);
			renderer->RemoveActor(actorEdge4);
			renderer->RemoveActor(actorEdge5);
			renderer->AddActor(actorEdge2);
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			break;
		case '3':
			renderer->RemoveActor(actorEdge1);
			renderer->RemoveActor(actorEdge2);
			renderer->RemoveActor(actorEdge3);
			renderer->RemoveActor(actorEdge4);
			renderer->RemoveActor(actorEdge5);
			renderer->AddActor(actorEdge3);
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			break;
		case '4':
			renderer->RemoveActor(actorEdge1);
			renderer->RemoveActor(actorEdge2);
			renderer->RemoveActor(actorEdge3);
			renderer->RemoveActor(actorEdge4);
			renderer->RemoveActor(actorEdge5);
			renderer->AddActor(actorEdge4);
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			break;
		case '5':
			renderer->RemoveActor(actorEdge1);
			renderer->RemoveActor(actorEdge2);
			renderer->RemoveActor(actorEdge3);
			renderer->RemoveActor(actorEdge4);
			renderer->RemoveActor(actorEdge5);
			renderer->AddActor(actorEdge5);
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			break;
		case '7':
			renderer->RemoveActor(actorPoly1);
			renderer->RemoveActor(actorPoly2);
			renderer->RemoveActor(actorPoly3);
			renderer->AddActor(actorPoly1);
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			break;
		case '8':
			renderer->RemoveActor(actorPoly1);
			renderer->RemoveActor(actorPoly2);
			renderer->RemoveActor(actorPoly3);
			renderer->AddActor(actorPoly2);
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			break;
		case '9':
			renderer->RemoveActor(actorPoly1);
			renderer->RemoveActor(actorPoly2);
			renderer->RemoveActor(actorPoly3);
			renderer->AddActor(actorPoly3);
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			break;
		case 't':
			{
			modTruncate.Step();
			originalTruncate.Step();
			vtkSmartPointer<vtkActor> edge_actor3 = CreateStrightUpPipeline(modTruncate.GetGraph());
			if (actorEdge3)
				actorEdge3->ShallowCopy(edge_actor3);
			else
				actorEdge3 = edge_actor3;
			vtkSmartPointer<vtkActor> edge_actor4 = CreateStrightUpPipeline(originalTruncate.GetGraph());
			if (actorEdge4)
				actorEdge4->ShallowCopy(edge_actor4);
			else
				actorEdge4 = edge_actor4;
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			}
			break;
		case 'g':
			{
			modTruncate.Process();
			originalTruncate.Process();
			cout << "num original edge:" << originalTruncate.GetGraph()->GetNumberOfEdges() << endl;
			vtkSmartPointer<vtkActor> edge_actor3 = CreateStrightUpPipeline(modTruncate.GetGraph());
			if (actorEdge3)
				actorEdge3->ShallowCopy(edge_actor3);
			else
				actorEdge3 = edge_actor3;
			vtkSmartPointer<vtkActor> edge_actor4 = CreateStrightUpPipeline(originalTruncate.GetGraph());
			if (actorEdge4)
				actorEdge4->ShallowCopy(edge_actor4);
			else
				actorEdge4 = edge_actor4;
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			double time_original = exact_algorithm->GetConsumedTime() + originalTruncate.GetTimeConsumed();
			double time_proposed = exact_algorithm->GetConsumedTime() + exact_algorithm->GetConsumedTime2() +  modTruncate.GetTimeConsumed();
			cout << "time consume original:" << time_original << "sec" << endl;
			cout << "time consume proposed:" << time_proposed << "sec" << endl;
			}
			break;
		case 'c':
			//cutting 
			if (modTruncate.isReadyToCut())
			{
				std::string filename;
				if (GetFileName(filename,"PLY File\0*.ply\0All Files\0*.*\0",true))
				{
					vtkSmartPointer<vtkPolyData> outputPropose = modTruncate.GetDiskTopologyPolydata();
					vtkSmartPointer<vtkPolyData> outputOriginal = originalTruncate.GetDiskTopologyPolydata();
					vtkSmartPointer<vtkPLYWriter> plyWriter = vtkSmartPointer<vtkPLYWriter>::New();
					
					std::string propose_filename = filename + std::string("_propose.ply");
					std::string original_filename = filename +  std::string("_original.ply");
					
					plyWriter->SetInputData(outputPropose);
					plyWriter->SetFileName(propose_filename.c_str());
					plyWriter->Update();
					cout << "Disk topology mesh export :" << propose_filename.c_str() << endl;

					plyWriter->SetInputData(outputOriginal);
					plyWriter->SetFileName(original_filename.c_str());
					plyWriter->Update();					
					cout << "Disk topology mesh export :" << original_filename.c_str() << endl;
				}
			}
			else
				cout << "not ready to export disk topoly mesh." << endl;
			break;
		
			break;
	}
		
	std::string key(ch);
	if (key == "F10")
	{
		ScreenShot(iren->GetRenderWindow());
	}
	else if (key == "F6")
	{
		//homotopy cutting		
		modTruncate.Process();
		originalTruncate.Process();
		cout << "num original edge:" << originalTruncate.GetGraph()->GetNumberOfEdges() << endl;
		vtkSmartPointer<vtkActor> edge_actor3 = CreateStrightUpPipeline(modTruncate.GetGraph());
		if (actorEdge3)
			actorEdge3->ShallowCopy(edge_actor3);
		else
			actorEdge3 = edge_actor3;
		vtkSmartPointer<vtkActor> edge_actor4 = CreateStrightUpPipeline(originalTruncate.GetGraph());
		if (actorEdge4)
			actorEdge4->ShallowCopy(edge_actor4);
		else
			actorEdge4 = edge_actor4;
		renderer->Modified();
		iren->GetRenderWindow()->Render();
		double time_original = exact_algorithm->GetConsumedTime() + originalTruncate.GetTimeConsumed();
		double time_proposed = exact_algorithm->GetConsumedTime() + exact_algorithm->GetConsumedTime2() +  modTruncate.GetTimeConsumed();
		cout << "Homotopy Cutting finished..."<< endl;
		cout << "time consume original:" << time_original << "sec" << endl;
		cout << "time consume proposed:" << time_proposed << "sec" << endl;
		cout << "========================================" << endl;
	}
	else if (key == "F7")
	{
		//iterated augment cutting
		double time_original = 0;
		if (originalTruncate.isReadyToCut())
		{
			cout << "Original Iterated Augment Cutting Started..."<< endl;
			vtkSmartPointer<vtkPolyData> outputOriginal = originalTruncate.GetDiskTopologyPolydata();
			CPolygonsData polygon ;
			polygon.InitailDiskTopology(outputOriginal);
			polygon.IteratedAugmentCutOriginal(&time_original,NULL);
			cout << "Original Iterated Augment Cutting Finished..."<< endl;
		}


		double time_proposed = 0;
		if (modTruncate.isReadyToCut())
		{
			cout << "Proposed Iterated Augment Cutting Started..."<< endl;
			vtkSmartPointer<vtkPolyData> outputPropose = modTruncate.GetDiskTopologyPolydata();
			CPolygonsData polygon ;
			vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
			polygon.InitailDiskTopology(outputPropose);
			polygon.IteratedAugmentCut(&time_proposed,output.GetPointer());

			cout << "Proposed Iterated Augment Cutting Finished..."<< endl;
			cout << "num cell :" << output->GetNumberOfPolys() << endl;


			

		}

		

		cout << "Iterated Augment Cuttings finished..."<< endl;
		cout << "time consume original:" << time_original << "sec" << endl;
		cout << "time consume proposed:" << time_proposed << "sec" << endl;
		cout << "========================================" << endl;
		
	}
	else if (key == "F8")
	{
		//Parameterization
		
		if (modTruncate.isReadyToCut())
		{
			vtkSmartPointer<vtkPolyData> outputPropose = modTruncate.GetDiskTopologyPolydata();
			CPolygonsData polygon ;
			polygon.InitailDiskTopology(outputPropose);
			polygon.Parameterize();
		}
			
	}
	else if (key == "plus")
	{
		double opac = actorMainPoly->GetProperty()->GetOpacity();
		opac = std::min(opac + 0.05, 1.0);
		actorMainPoly->GetProperty()->SetOpacity(opac);
		renderer->Modified();
		iren->GetRenderWindow()->Render();
	}
	else if (key == "minus")
	{
		double opac = actorMainPoly->GetProperty()->GetOpacity();
		opac = std::max(opac - 0.05, 0.05);
		actorMainPoly->GetProperty()->SetOpacity(opac);
		renderer->Modified();
		iren->GetRenderWindow()->Render();
	}
	//std::cout << "Pressed: " << iren->GetKeySym() << endl;
}
void pickCallbackFunc(vtkObject* caller, unsigned long eid, void* clientdata, void *calldata)
{
	
	vtkPointPicker* picker = static_cast<vtkPointPicker*>(caller);
	MouseInteractorStylePP* TrackballStyle = static_cast<MouseInteractorStylePP*>(clientdata);

	if (picker->GetPointId() >= 0 )
	{
		int vertexID = picker->GetPointId();
		std::cout << "Pick Point:" << vertexID <<  std::endl;

		if (TrackballStyle->GetInteractor()->GetControlKey() != 0)
		{
			vtkSmartPointer<vtkPolyData> PolyData = polydata;
			Process(PolyData,picker->GetPointId());
			ColoredPoint( TrackballStyle->GetInteractor()->GetRenderWindow()->GetRenderers()->GetFirstRenderer(),PolyData->GetPoint(vertexID), 0.0,1.0,0.0);
			actorEdge1->Modified();
			actorEdge2->Modified();
			TrackballStyle->GetInteractor()->GetRenderWindow()->Render();
		}
	}
	

}
void LoadToGeodesic(vtkSmartPointer<vtkPolyData> polydata)
{
	const int numPoints = polydata->GetNumberOfPoints();
	const int numFaces = polydata->GetPolys()->GetNumberOfCells();
	
	std::vector<double> vertices(numPoints*3);
	std::vector<int> faces(numFaces*3);

	for (vtkIdType i = 0 ; i < numPoints; i++)
	{
		polydata->GetPoint(i,&vertices[i*3]);		
	}
	std::vector<int>::iterator faceItr = faces.begin();
	vtkSmartPointer<vtkCellArray> cellArray = polydata->GetPolys();
	cellArray->InitTraversal();
	vtkIdType numPts;
	vtkIdType *pts;
	while (cellArray->GetNextCell(numPts,pts))
	{
		*faceItr++ = (int)pts[0];
		*faceItr++ = (int)pts[1];
		*faceItr++ = (int)pts[2];
	}
	
	geosesic_mesh.clear_memory();
	geosesic_mesh.initialize_mesh_data(vertices,faces);
}


bool GenerateGeodesicDistance(geodesic::GeodesicAlgorithmExact &geoAlgo,int vertexID , std::vector<collisionEdgeInfo>& collision_edges)
{
	
	if (vertexID < 0 || vertexID >= geoAlgo.mesh()->vertices().size())
		return false;
	
	geodesic::SurfacePoint source(&geoAlgo.mesh()->vertices()[vertexID]);		//create source 
	std::vector<geodesic::SurfacePoint> all_sources(1,source);					//in general, there could be multiple sources, but now we have only one
	
	geoAlgo.propagate(all_sources);	//cover the whole mesh

	
	collision_edges.clear();
	for (int i = 0 ; i < geoAlgo.CollisionEdges().size() ; i++)
	{
		geodesic::edge_pointer e = geoAlgo.CollisionEdges()[i];	
		collisionEdgeInfo info;
		info.v0 = e->adjacent_vertices()[0]->id();
		info.v1 = e->adjacent_vertices()[1]->id();
		info.type = e->source_collision();
		collision_edges.push_back(info);
	}
	return true;
}

vtkSmartPointer<vtkMutableUndirectedGraph> GenerateCollisionEdgeGraph(vtkSmartPointer<vtkPolyData> source_polydata, std::vector<collisionEdgeInfo>& collision_edges)
{
	vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	for(vtkIdType i = 0; i < source_polydata->GetNumberOfPoints(); i++)
    {
		graph->AddVertex();
    }
  
	vtkSmartPointer<vtkPoints> points = 
    vtkSmartPointer<vtkPoints>::New();
	points->ShallowCopy(source_polydata->GetPoints());
	/*
	for (int p = 0 ; p < source_polydata->GetNumberOfPoints(); p++)
		points->InsertNextPoint(source_polydata->GetPoint(p));
	*/
 
	// Add the coordinates of the points to the graph
	graph->SetPoints(points);
 
	//graph->GetVertexData()->PassData(source_polydata->GetPointData());
	int numCollisionEdge = (int) collision_edges.size();
	for (int i = 0 ; i < numCollisionEdge ; i++)
	{
		if (collision_edges[i].type == 1|| collision_edges[i].type == 2)
			graph->AddEdge(collision_edges[i].v0,collision_edges[i].v1);
	}
	return graph;
}

vtkSmartPointer<vtkMutableUndirectedGraph> TruncateGraph(vtkSmartPointer<vtkMutableUndirectedGraph> source_graph)
{
	vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	graph->DeepCopy(source_graph);
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
				/*
				vtkSmartPointer<vtkAdjacentVertexIterator> adjacent =vtkSmartPointer<vtkAdjacentVertexIterator>::New();
				graph->GetAdjacentVertices(vid[i],adjacent);
				int adjnum = 0;
				while (adjacent->HasNext())
				{
					vtkIdType u = adjacent->Next();
					++adjnum;
				}
				*/
				
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
	while (numRemoveEdges > 0);
	int numEdges = graph->GetNumberOfEdges();
	return graph;
}

vtkSmartPointer<vtkMutableUndirectedGraph> StraightupGraph(vtkSmartPointer<vtkMutableUndirectedGraph> source_graph, OmMesh& mesh)
{
	vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	graph->DeepCopy(source_graph);

	int numRemoveEdges = 0;

	do
	{
		numRemoveEdges = 0;
		vtkSmartPointer<vtkIdTypeArray> removeEdgesArray = vtkSmartPointer<vtkIdTypeArray>::New();
		int numEdges = graph->GetNumberOfEdges();
		vtkSmartPointer<vtkEdgeListIterator> eit =vtkSmartPointer<vtkEdgeListIterator>::New();
		graph->GetEdges(eit);

		while (eit->HasNext())
		{
			vtkEdgeType e = eit->Next();
			vtkIdType vid[2] = {e.Source,e.Target};
		
			for (int i = 0 ; i < 2 ;i++)
			{
				int adjnum = graph->GetDegree(vid[i]);
				vtkIdType adjvid[2];
				if (adjnum > 2)
				{
					continue;
				}
				else
				{
					vtkSmartPointer<vtkAdjacentVertexIterator> adjacent =vtkSmartPointer<vtkAdjacentVertexIterator>::New();
					graph->GetAdjacentVertices(vid[i],adjacent);
					adjvid[0] = adjacent->Next();
					adjvid[1] = adjacent->Next();
				}
			
				//check for mutual face share 
			
				OmMesh::HalfedgeHandle heh;
				heh = mesh.find_halfedge(mesh.vertex_handle(adjvid[0]),mesh.vertex_handle(adjvid[1]));

				if (!heh.is_valid())
				{
					heh = mesh.find_halfedge(mesh.vertex_handle(adjvid[1]),mesh.vertex_handle(adjvid[0]));
				}
				

				if (heh.is_valid())
				{
					//share face!
					if (graph->GetEdgeId(adjvid[0],adjvid[1]) < 0)
					{
						graph->AddEdge(adjvid[0],adjvid[1]);
						vtkIdType eid0 = graph->GetEdgeId(vid[i],adjvid[0]);
						vtkIdType eid1 = graph->GetEdgeId(vid[i],adjvid[1]);					
						removeEdgesArray->InsertNextValue(eid0);					
						removeEdgesArray->InsertNextValue(eid1);
						graph->RemoveEdges(removeEdgesArray);
						numRemoveEdges = 2;
						break;
					}
				}
			}

			if (numRemoveEdges > 0)
				break;
		}			
	}
	while (numRemoveEdges > 0);
	vtkSmartPointer<vtkMutableUndirectedGraph> newGraph = TruncateGraph(graph);
	graph->ShallowCopy(newGraph);
	return graph;


}

vtkSmartPointer<vtkMutableUndirectedGraph> CleanGraph(vtkSmartPointer<vtkMutableUndirectedGraph> source_graph, OmMesh &mesh )
{
	vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	graph->DeepCopy(source_graph);
	int numEdgesBefore = 0;
	
	//clean parallel edge paths	
	do
	{
		numEdgesBefore = graph->GetNumberOfEdges();
		vtkSmartPointer<vtkEdgeListIterator> eit =vtkSmartPointer<vtkEdgeListIterator>::New();
		graph->GetEdges(eit);
		
		while (eit->HasNext())
		{
			vtkEdgeType e = eit->Next();
			vtkIdType vid[2*3] = {e.Source,e.Target};
			if (vid[0] > vid[1])
			{
				vid[0] = e.Target;
				vid[1] = e.Source;
			}
			
			bool val2 = true;
			for (int i = 0 ; i < 2 ; i++)
			{
				int adjnum = graph->GetDegree(vid[i]);
				if (adjnum > 2)
				{
					//not all valence-2 vertex edge
					val2 = false;
					break;
				}
			} 			
			if (!val2)
				continue;
			
			
			for (int i = 0 ; i< 2; i++)
			{
				vtkSmartPointer<vtkAdjacentVertexIterator> adjacent =vtkSmartPointer<vtkAdjacentVertexIterator>::New();
				graph->GetAdjacentVertices(vid[i],adjacent);
				vid[2 + i*2 + 0] = vid[i];
				
				
				vid[2 + i*2 + 1] = adjacent->Next();
				if (vid[2 + i*2 + 1] == vid[(i+1)%2])
					vid[2 + i*2 + 1] = adjacent->Next();
			}

			bool delete_edge = true;

			for (int i = 0 ; i < 3 ;i++) //3 edges
			{
				//check if radius of opposite vertex has vertex in cutting path
				OmMesh::VertexHandle vh[2] = {mesh.vertex_handle(vid[i*2 + 0]), mesh.vertex_handle(vid[i*2 + 1])};
				OmMesh::HalfedgeHandle heh[2] ;
				heh[0] = mesh.find_halfedge(vh[0],vh[1]);
				if (!heh[0].is_valid())
				{
					heh[0] = mesh.find_halfedge(vh[1],vh[0]);
				}
				else
				{
					heh[1] = mesh.opposite_halfedge_handle(heh[0]);
				}
			
				bool hitted = false;
				for (int i = 0; i < 2; i ++)
				{
					if (!heh[i].is_valid())
						continue;
					int faceid = mesh.face_handle(heh[i]).idx(); 
					int vertexid = mesh.opposite_vh(heh[i]).idx();

					//search vertexid in graph
					if (graph->GetDegree(vertexid) >= 2)
					{
						hitted = true;					
						break;
					}				
				}
				if (!hitted)
				{
					delete_edge = false;
					break;
				}
			}
			
			if (delete_edge)
			{
				graph->RemoveEdge(e.Id);
				vtkSmartPointer<vtkMutableUndirectedGraph> newGraph = TruncateGraph(graph);
				newGraph = StraightupGraph(newGraph,mesh);
				graph->ShallowCopy(newGraph);
				break;					
			}
		}
	}
	while (numEdgesBefore > graph->GetNumberOfEdges());

	//return graph;
	//clean small cycle edge paths
	do
	{
		numEdgesBefore = graph->GetNumberOfEdges();
		vtkSmartPointer<vtkEdgeListIterator> eit =vtkSmartPointer<vtkEdgeListIterator>::New();
		graph->GetEdges(eit);
		
		while (eit->HasNext())
		{
			bool remove_edge0 = false;
			vtkEdgeType e = eit->Next();
			vtkIdType vid[2] = {e.Source,e.Target};
			bool val2 = true;
			for (int i = 0 ; i < 2 ; i++)
			{
				int adjnum = graph->GetDegree(vid[i]);
				if (adjnum >= 3)
				{

					OmMesh::VertexHandle v0 = mesh.vertex_handle(vid[i]);
					OmMesh::HalfedgeHandle heh0 = mesh.halfedge_handle(v0);
					OmMesh::HalfedgeHandle heh1;
					OmMesh::HalfedgeHandle startheh = heh0;
					
					vtkIdType vtkEdgeID0,vtkEdgeID1,vtkEdgeID2,vtkEdgeID3,vtkEdgeID4= 0;
					/*
					//        v2 
					//      /    \
					//    -v0  | v4-
					//      \    /
					//        v1
					*/
					do 
					{
						OmMesh::VertexHandle v1 = mesh.to_vertex_handle(heh0);
						heh1 = mesh.next_halfedge_handle(heh0);
						OmMesh::VertexHandle v2 = mesh.to_vertex_handle(heh1);
						if ((vtkEdgeID0=graph->GetEdgeId(v0.idx(),v1.idx())) > 0)
						{
							if ((vtkEdgeID1=graph->GetEdgeId(v0.idx(),v2.idx())) > 0)
							{
								
								if ((vtkEdgeID2=graph->GetEdgeId(v1.idx(),v2.idx())) > 0)
								{
									//triangle cycle
									//remove heh0
									remove_edge0 = true;
								}
								else
								{
									OmMesh::HalfedgeHandle opheh = mesh.opposite_halfedge_handle(heh1);
									if (opheh.is_valid())
									{
										OmMesh::VertexHandle v4 = mesh.opposite_vh(opheh);
										if (((vtkEdgeID3=graph->GetEdgeId(v4.idx(),v1.idx())) > 0) &&
											((vtkEdgeID4=graph->GetEdgeId(v4.idx(),v2.idx())) > 0))
										{
											//rectangle cycle
											//remove heh0
											remove_edge0 = true;
										}
									}
								}

								if (remove_edge0)
								{
									graph->RemoveEdge(vtkEdgeID0);
									vtkSmartPointer<vtkMutableUndirectedGraph> newGraph = TruncateGraph(graph);
									
									newGraph = StraightupGraph(newGraph,mesh);
									graph->ShallowCopy(newGraph);
									
									break;
								}
							}
						}
						heh0 =  mesh.opposite_halfedge_handle(mesh.next_halfedge_handle(heh1));
					}
					while (heh0.is_valid() && heh0 != startheh);
					
				}
				if (remove_edge0)
					break;
			}
			if (remove_edge0)
				break;
		}
	}
	while (numEdgesBefore > graph->GetNumberOfEdges());
	return graph;
}

vtkSmartPointer<vtkActor> CreateBeforeTruncatePipeline(vtkSmartPointer<vtkMutableUndirectedGraph> BTGraph,std::vector<collisionEdgeInfo>& collision_edges)
{
	vtkSmartPointer<vtkGraphToPolyData> g2pBT = vtkSmartPointer<vtkGraphToPolyData>::New(); 
	g2pBT->SetInputData(BTGraph);
	g2pBT->Update();
	vtkSmartPointer<vtkPolyData> collisionEdgesBTPolydata = vtkSmartPointer<vtkPolyData>::New();
	collisionEdgesBTPolydata->ShallowCopy(g2pBT->GetOutput());
	
	


	// Setup the colors array
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(4);
	colors->SetName("Colors");
	// Setup two colors 
	unsigned char white[4] = {255, 255, 255,0};
	unsigned char red[4] = {255, 0, 0,0};
	unsigned char green[4] = {0, 255, 0,255};
	unsigned char blue[4] = {0, 0, 255,255};
	unsigned char yellow[4] = {255, 255, 0,255};
	int numEdge = collisionEdgesBTPolydata->GetLines()->GetNumberOfCells();
	for (int cellID = 0 ; cellID < collisionEdgesBTPolydata->GetLines()->GetNumberOfCells(); cellID++)
	{
		if (collision_edges[cellID].type == 1)
			colors->InsertNextTupleValue(blue);	
		else if (collision_edges[cellID].type == 2)
			colors->InsertNextTupleValue(green);

	}

	

	/*
	//special edges former FROM_BOTHFACE tag
	int numCollisionEdge = (int) collision_edges.size();
	for (int i = 0 ; i < numCollisionEdge ; i++)
	{
		if (collision_edges[i].type == -1)
		{
			vtkIdType eid[2] = {collision_edges[i].v0,collision_edges[i].v1};			

			vtkIdType newID = collisionEdgesBTPolydata->GetLines()->InsertNextCell(2,eid);
			
			colors->InsertNextTupleValue(yellow);
		}
	}
	*/
	collisionEdgesBTPolydata->GetCellData()->SetScalars(colors);

	vtkSmartPointer<vtkPolyDataMapper> edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	edge_mapper->SetInputData(collisionEdgesBTPolydata);
 
	vtkSmartPointer<vtkActor> edge_actor = vtkSmartPointer<vtkActor>::New();
	edge_actor->SetMapper(edge_mapper);
	edge_actor->GetProperty()->SetLineWidth(5.0);
	edge_actor->GetProperty()->SetOpacity(1.0);
	return edge_actor;
}
vtkSmartPointer<vtkActor> CreateAfterTruncatePipeline(vtkSmartPointer<vtkMutableUndirectedGraph> ATGraph,std::vector<collisionEdgeInfo>& collision_edges)
{
	vtkSmartPointer<vtkGraphToPolyData> g2pAT = vtkSmartPointer<vtkGraphToPolyData>::New(); 
	g2pAT->SetInputData(ATGraph);
	g2pAT->Update();
	vtkSmartPointer<vtkPolyData> collisionEdgesATPolydata = vtkSmartPointer<vtkPolyData>::New();
	collisionEdgesATPolydata->ShallowCopy(g2pAT->GetOutput());

	// Setup the colors array
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");
	// Setup two colors 
	unsigned char red[3] = {255, 0, 0};
	unsigned char green[3] = {0, 255, 0};
	unsigned char blue[3] = {0, 0, 255};
	unsigned char yellow[3] = {255, 255, 0};
	int numEdge = collisionEdgesATPolydata->GetLines()->GetNumberOfCells();
	for (int cellID = 0 ; cellID < collisionEdgesATPolydata->GetLines()->GetNumberOfCells(); cellID++)
	{
		if (collision_edges[cellID].type == 1)
			colors->InsertNextTupleValue(blue);	
		else if (collision_edges[cellID].type == 2)
			colors->InsertNextTupleValue(green);
		
	}
	collisionEdgesATPolydata->GetCellData()->SetScalars(colors);


	vtkSmartPointer<vtkPolyDataMapper> edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	edge_mapper->SetInputData(collisionEdgesATPolydata);
 
	vtkSmartPointer<vtkActor> edge_actor = vtkSmartPointer<vtkActor>::New();
	edge_actor->SetMapper(edge_mapper);
	edge_actor->GetProperty()->SetLineWidth(4.0);
	edge_actor->GetProperty()->SetOpacity(1.0);
	return edge_actor;
}

vtkSmartPointer<vtkActor> CreateStrightUpPipeline(vtkSmartPointer<vtkMutableUndirectedGraph> SUGraph)
{
	vtkSmartPointer<vtkMutableUndirectedGraph> modGraph =  vtkSmartPointer<vtkMutableUndirectedGraph>::New();


	vtkSmartPointer<vtkGraphToPolyData> g2pAT = vtkSmartPointer<vtkGraphToPolyData>::New(); 
	g2pAT->SetInputData(SUGraph);
	g2pAT->Update();
	vtkSmartPointer<vtkPolyData> suPolydata = vtkSmartPointer<vtkPolyData>::New();
	suPolydata->ShallowCopy(g2pAT->GetOutput());
	
	vtkSmartPointer<vtkUnsignedCharArray> edgecolors = vtkUnsignedCharArray::SafeDownCast( SUGraph->GetEdgeData()->GetArray("EdgeColors"));
	vtkSmartPointer<vtkUnsignedCharArray> collisionTags = vtkUnsignedCharArray::SafeDownCast( SUGraph->GetEdgeData()->GetArray("CollisionTag"));

	// Setup the colors array
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(4);
	colors->SetName("Colors");
	// Setup two colors 
	unsigned char red[3] = {255, 0, 0};
	unsigned char green[3] = {0, 255, 0};
	unsigned char blue[3] = {0, 0, 255};
	unsigned char yellow[3] = {255, 255, 0};
	int numEdge = suPolydata->GetLines()->GetNumberOfCells();
	for (int cellID = 0 ; cellID < suPolydata->GetLines()->GetNumberOfCells(); cellID++)
	{
		
		unsigned char edgeColor[4];


		edgecolors->GetTupleValue(cellID,edgeColor);
		colors->InsertNextTupleValue(edgeColor);
		
	}
	suPolydata->GetCellData()->SetScalars(colors);
	

	vtkSmartPointer<vtkPolyDataMapper> edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	edge_mapper->SetInputData(suPolydata);
 
	vtkSmartPointer<vtkActor> edge_actor = vtkSmartPointer<vtkActor>::New();
	edge_actor->SetMapper(edge_mapper);
	edge_actor->GetProperty()->SetLineWidth(4.0);
	edge_actor->GetProperty()->SetOpacity(1.0);
	return edge_actor;
}

vtkSmartPointer<vtkActor> CreateCleanPipeline(vtkSmartPointer<vtkMutableUndirectedGraph> CleanGraph)
{
	vtkSmartPointer<vtkGraphToPolyData> g2pAT = vtkSmartPointer<vtkGraphToPolyData>::New(); 
	g2pAT->SetInputData(CleanGraph);
	g2pAT->Update();
	vtkSmartPointer<vtkPolyData> cleanPolydata = vtkSmartPointer<vtkPolyData>::New();
	cleanPolydata->ShallowCopy(g2pAT->GetOutput());

	// Setup the colors array
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");
	// Setup two colors 
	unsigned char purple[3] = {163, 73, 164};
	int numEdge = cleanPolydata->GetLines()->GetNumberOfCells();
	for (int cellID = 0 ; cellID < cleanPolydata->GetLines()->GetNumberOfCells(); cellID++)
	{
		
		colors->InsertNextTupleValue(purple);
		
	}
	cleanPolydata->GetCellData()->SetScalars(colors);


	vtkSmartPointer<vtkPolyDataMapper> edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	edge_mapper->SetInputData(cleanPolydata);
 
	vtkSmartPointer<vtkActor> edge_actor = vtkSmartPointer<vtkActor>::New();
	edge_actor->SetMapper(edge_mapper);
	edge_actor->GetProperty()->SetLineWidth(1.0);
	edge_actor->GetProperty()->SetOpacity(0.75);
	return edge_actor;
}

vtkSmartPointer<vtkActor> CreateFinalPipeline(vtkSmartPointer<vtkMutableUndirectedGraph> FinalGraph)
{
	vtkSmartPointer<vtkGraphToPolyData> g2pAT = vtkSmartPointer<vtkGraphToPolyData>::New(); 
	g2pAT->SetInputData(FinalGraph);
	g2pAT->Update();
	vtkSmartPointer<vtkPolyData> finalPolydata = vtkSmartPointer<vtkPolyData>::New();
	finalPolydata->ShallowCopy(g2pAT->GetOutput());

	// Setup the colors array
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");
	// Setup two colors 
	unsigned char purple[3] = {163, 73, 164};
	int numEdge = finalPolydata->GetLines()->GetNumberOfCells();
	for (int cellID = 0 ; cellID < finalPolydata->GetLines()->GetNumberOfCells(); cellID++)
	{
		
		colors->InsertNextTupleValue(purple);
		
	}
	finalPolydata->GetCellData()->SetScalars(colors);


	vtkSmartPointer<vtkPolyDataMapper> edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	edge_mapper->SetInputData(finalPolydata);
 
	vtkSmartPointer<vtkActor> edge_actor = vtkSmartPointer<vtkActor>::New();
	edge_actor->SetMapper(edge_mapper);
	edge_actor->GetProperty()->SetLineWidth(2.0);
	edge_actor->GetProperty()->SetOpacity(0.75);
	return edge_actor;
}

void ColoredPoint(vtkSmartPointer<vtkRenderer> _renderer ,  double pt[3],double r, double g, double b)
{
	vtkSmartPointer<vtkPoints> point = vtkSmartPointer<vtkPoints>::New();
	point->InsertNextPoint(pt);
	
	if (source_point)
	{
		source_point->ShallowCopy(point);
		source_point->Modified();
	}
	else
	{
		source_point = point;
		vtkSmartPointer<vtkPolyData> pointsPolydata =  vtkSmartPointer<vtkPolyData>::New();
		pointsPolydata->SetPoints(source_point);
		vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =  vtkSmartPointer<vtkVertexGlyphFilter>::New();
		vertexGlyphFilter->AddInputData(pointsPolydata);
		vertexGlyphFilter->Update();
 
		vtkSmartPointer<vtkPolyDataMapper> pointsMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		pointsMapper->SetInputConnection(vertexGlyphFilter->GetOutputPort());

		vtkSmartPointer<vtkActor> pointsActor =  vtkSmartPointer<vtkActor>::New();
		pointsActor->SetMapper(pointsMapper);
		pointsActor->GetProperty()->SetPointSize(15);
		pointsActor->GetProperty()->SetColor(r,g,b);

		_renderer->AddActor(pointsActor);
	}
}

void GetNeighborCell(vtkSmartPointer<vtkPolyData> source_polydata, vtkIdType vid , int level, std::vector<vtkIdType> &storeIDs)
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

vtkSmartPointer<vtkPolyData> CreateSurroundGraphPolydata(vtkSmartPointer<vtkMutableUndirectedGraph> source_graph, vtkSmartPointer<vtkPolyData> source_polydata,int level)
{
	if (level < 1)
		return NULL;
	vtkSmartPointer<vtkEdgeListIterator> eit =vtkSmartPointer<vtkEdgeListIterator>::New();
	source_graph->GetEdges(eit);
	vtkSmartPointer<vtkPolyData> outputPoly = vtkSmartPointer<vtkPolyData>::New();
	outputPoly->DeepCopy(source_polydata);
	source_polydata->BuildLinks();
	bool *cellCheck = new bool[source_polydata->GetPolys()->GetNumberOfCells()];
	memset(cellCheck,0,sizeof(bool)*source_polydata->GetPolys()->GetNumberOfCells());

	

	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
	while (eit->HasNext())
	{
		vtkEdgeType e = eit->Next();
		vtkIdType vid[2] = {e.Source,e.Target};
		
		for (int i = 0 ; i < 2 ; i++)
		{
			std::vector<vtkIdType> add_cellIDs;
			GetNeighborCell(source_polydata,vid[i],level,add_cellIDs);
			for (int j = 0 ; j < add_cellIDs.size(); j++)
			{
				if (!cellCheck[add_cellIDs[j]])
				{
					vtkCell *cell = source_polydata->GetCell(add_cellIDs[j]);
				
					cellArray->InsertNextCell(cell);
					cellCheck[add_cellIDs[j]] = true;
				}
			}
		}		
	}

	outputPoly->SetPolys(cellArray);
	outputPoly->BuildLinks();
	delete [] cellCheck;
	return outputPoly;

}



bool seedRemoved = false;
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
vtkSmartPointer<vtkMutableUndirectedGraph> GIMtruncate(vtkSmartPointer<vtkPolyData> polydata, geodesic::GeodesicAlgorithmExact &geo)
{
	OmMesh mesh;
	vtkPolydata2OpenMesh(polydata,&mesh );	

	vtkSmartPointer<vtkMutableUndirectedGraph> graph = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	for(vtkIdType i = 0; i < polydata->GetNumberOfPoints(); i++)
    {
		graph->AddVertex();
    }
  
	vtkSmartPointer<vtkPoints> points = 
    vtkSmartPointer<vtkPoints>::New();
	points->ShallowCopy(polydata->GetPoints());
	graph->SetPoints(points);
	for (OmMesh::EdgeIter e_it=mesh.edges_begin(); e_it!=mesh.edges_end(); ++e_it) 
	{
		OmMesh::EdgeHandle eh = *e_it;
		if (!mesh.is_boundary(eh))
		{
			OmMesh::HalfedgeHandle heh = mesh.halfedge_handle(eh,0);
			if (!heh.is_valid())
				heh = mesh.halfedge_handle(eh,1);
			OmMesh::HalfedgeHandle heh_prev = mesh.prev_halfedge_handle(heh);
			int a = mesh.to_vertex_handle(heh).idx();
			int b = mesh.to_vertex_handle(heh_prev).idx();
			graph->AddEdge(a,b);		
		}
	}

	//return graph;
	//remove seed triangle (nearest (geodesic) from source point)
	
	OmMesh::FaceHandle seedFace;
	double distanceSeed = geodesic::GEODESIC_INF;
	
	for (OmMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end(); ++f_it)
	{
		OmMesh::FaceHandle fh = *f_it;
		OmMesh::HalfedgeHandle heh = mesh.halfedge_handle(fh);
		unsigned int pts[3];
		pts[0] = mesh.to_vertex_handle(heh).idx();
		
		heh = mesh.next_halfedge_handle(heh);
		pts[1] = mesh.to_vertex_handle(heh).idx();

		heh = mesh.next_halfedge_handle(heh);
		pts[2] = mesh.to_vertex_handle(heh).idx();

		geodesic::face_pointer faceG =  geo.mesh()->find_face(pts[0],pts[1],pts[2]);
		if (faceG)
		{
			double best_distance;
			geodesic::SurfacePoint p(faceG);
			geo.best_source(p,best_distance);
			if (distanceSeed > best_distance)
			{
				seedFace = fh;
				distanceSeed = best_distance;
				
			}		
		}
	}

	std::map<double,hedge_data> candidate_edges;
	//remove seed triangle
	{
		OmMesh::FaceHandle removeFaceHandle = seedFace;
		OmMesh::HalfedgeHandle candidate_heh[3];
			
		candidate_heh[0] = mesh.halfedge_handle(removeFaceHandle);
		candidate_heh[1] = mesh.next_halfedge_handle(candidate_heh[0]);
		candidate_heh[2] = mesh.next_halfedge_handle(candidate_heh[1]);
		for (int i = 0 ; i < 3 ; i++)
		{
			vtkIdType endID,startID ;
			endID = mesh.to_vertex_handle(candidate_heh[i]).idx();
			startID = mesh.to_vertex_handle( mesh.prev_halfedge_handle(candidate_heh[i])).idx();
			vtkIdType edgeID =  graph->GetEdgeId(endID,startID);
			
			if (!mesh.is_boundary(mesh.edge_handle(candidate_heh[i])))
			{		
				geodesic::edge_pointer edgeG = geo.mesh()->find_edge(endID,startID);
				double best_distance;
				geodesic::SurfacePoint p(edgeG);
				geo.best_source(p,best_distance);	
				hedge_data hedata(endID,startID); //reserve input because opposite hedge
				candidate_edges.insert(std::pair<double,hedge_data>(best_distance,hedata));
			}
			else
			{
				//graph->RemoveEdge(edgeID);
			}		
		}
		polydata->DeleteCell(seedFace.idx());
		mesh.delete_face(seedFace,false);
		mesh.garbage_collection();
		polydata->RemoveDeletedCells();
	}

	return graph;

	
	//remove non-boundary edge that adjacent only triangle
	int count = 0;
	while (candidate_edges.size() > 0 && count < 6)
	{
		count ++;
		hedge_data checkingEdge = candidate_edges.begin()->second;
		OmMesh::HalfedgeHandle checkingHeh = mesh.find_halfedge( mesh.vertex_handle(checkingEdge.startID) , mesh.vertex_handle(checkingEdge.endID));
		if (!checkingHeh.is_valid())
		{
			//this edge (of graph) have adjacent 0 edge
			//remove from map -> will be cutpath
			int aaa = 0;
			
		}
		else
		{
			//remove face and edge
			OmMesh::FaceHandle removeFaceHandle = mesh.face_handle(checkingHeh);
			OmMesh::HalfedgeHandle candidate_heh[2];
			
			candidate_heh[0] = mesh.next_halfedge_handle(checkingHeh);
			candidate_heh[1] = mesh.next_halfedge_handle(candidate_heh[0]);
			for (int i = 0 ; i < 2 ; i++)
			{
				vtkIdType endID,startID ;
				endID = mesh.to_vertex_handle(candidate_heh[i]).idx();
				startID = mesh.to_vertex_handle( mesh.prev_halfedge_handle(candidate_heh[i])).idx();
				vtkIdType edgeID =  graph->GetEdgeId(endID,startID);
				if (!mesh.is_boundary(mesh.edge_handle(candidate_heh[i])))
				{		
					geodesic::edge_pointer edgeG = geo.mesh()->find_edge(startID,endID);
					double best_distance;
					geodesic::SurfacePoint p(edgeG);
					geo.best_source(p,best_distance);		
					hedge_data hedata(endID,startID); //reserve input because of opposite hedge
					candidate_edges.insert(std::pair<double,hedge_data>(best_distance,hedata ));
				}
				else
				{
					//graph->RemoveEdge(edgeID);
				}		
			}

			vtkIdType edgeID =  graph->GetEdgeId(checkingEdge.startID,checkingEdge.endID);
			graph->RemoveEdge(edgeID);
			polydata->DeleteCell(removeFaceHandle.idx());
			mesh.delete_face( removeFaceHandle,false);
			mesh.garbage_collection();
			polydata->RemoveDeletedCells();
		}
		
		candidate_edges.erase( candidate_edges.begin());


	}


	return graph;
}

void ScreenShot(vtkRenderWindow* renderWindow)
{
	// Screenshot  

	std::string filename;
	if (GetFileName(filename,"PNG file\0*.png\0All\0*.*\0",true))
	{
		vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = 
		vtkSmartPointer<vtkWindowToImageFilter>::New();
		  windowToImageFilter->SetInput(renderWindow);
		  //windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
		  windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
		  windowToImageFilter->Update();
 
		  vtkSmartPointer<vtkPNGWriter> writer = 
			vtkSmartPointer<vtkPNGWriter>::New();
		  writer->SetFileName(filename.c_str());
		  writer->SetInputConnection(windowToImageFilter->GetOutputPort());
		  writer->Write();
		  
		  cout << "Screenshot :" << filename.c_str() << endl;
	}
  
}