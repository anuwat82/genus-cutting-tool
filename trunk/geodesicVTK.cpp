
#include "stdafx.h"

#include <Windows.h>
#include "utils.h"
#include "VTK_Header.h"
#include "MouseInteractorStylePP.h"
#include "vtkOFFReader.h"
#include <ostream>
#include <string>     // std::string, std::stoi
#include "geodesic\geodesic_algorithm_exact.h"
vtkWeakPointer<vtkPolyData> polydata;
vtkWeakPointer<vtkPoints> source_point;
	
vtkWeakPointer<vtkActor> actorEdge1;
vtkWeakPointer<vtkActor> actorEdge2;
geodesic::Mesh geosesic_mesh;

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
bool GenerateGeodesicDistance(int vertexID , std::vector<collisionEdgeInfo>& collision_edges);
vtkSmartPointer<vtkMutableUndirectedGraph> GenerateCollisionEdgeGraph(vtkSmartPointer<vtkPolyData> source_polydata, std::vector<collisionEdgeInfo>& collision_edges);
vtkSmartPointer<vtkMutableUndirectedGraph> TruncateGraph(vtkSmartPointer<vtkMutableUndirectedGraph> source_graph);

vtkSmartPointer<vtkActor> CreateBeforeTruncatePipeline(vtkSmartPointer<vtkMutableUndirectedGraph> BTGraph,std::vector<collisionEdgeInfo>& collision_edges);
vtkSmartPointer<vtkActor> CreateAfterTruncatePipeline(vtkSmartPointer<vtkMutableUndirectedGraph> ATGraph,std::vector<collisionEdgeInfo>& collision_edges);

void Process(vtkSmartPointer<vtkPolyData> polydata ,int sourceVertexID);
std::string GetFileExtension(const std::string& FileName)
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}

void Process(vtkSmartPointer<vtkPolyData> polydata , int sourceVertexID)
{
	std::vector<collisionEdgeInfo> collision_edges;
	
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


	GenerateGeodesicDistance(sourceVertexID,collision_edges);
	vtkSmartPointer<vtkMutableUndirectedGraph> collisionEdgesGraphBeforeTruncate = GenerateCollisionEdgeGraph(polydata,collision_edges);
	vtkSmartPointer<vtkMutableUndirectedGraph> collisionEdgesGraphAfterTruncate =  TruncateGraph(collisionEdgesGraphBeforeTruncate);
	
 
	vtkSmartPointer<vtkActor> edge_actor1 = CreateBeforeTruncatePipeline(collisionEdgesGraphBeforeTruncate, collision_edges); 
	vtkSmartPointer<vtkActor> edge_actor2 = CreateAfterTruncatePipeline(collisionEdgesGraphAfterTruncate, collision_edges);
	
	if (edge_actor1->GetReferenceCount() == 1)
		edge_actor1->SetReferenceCount(2);
	if (edge_actor2->GetReferenceCount() == 1)
		edge_actor2->SetReferenceCount(2);
	
	if (actorEdge1)
		actorEdge1->ShallowCopy(edge_actor1);
	else
		actorEdge1 = edge_actor1;
	if (actorEdge2)
		actorEdge2->ShallowCopy(edge_actor2);
	else
		actorEdge2 = edge_actor2;

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
		if (!GetModelFileName(filename))
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

	/*
	std::vector<collisionEdgeInfo> collision_edges;
	int numVertex = PLYReader->GetOutput()->GetNumberOfPoints();
	if (sourceVertex >= numVertex)
	{
		std::cout << "invalid source vertex id ... reset to id 0." << endl;
		sourceVertex = 0;
	}
	std::cout << "source vertex id " << sourceVertex << endl;

	LoadToGeodesic(PLYReader->GetOutput());
	GenerateGeodesicDistance(sourceVertex,collision_edges);
	
	vtkSmartPointer<vtkMutableUndirectedGraph> collisionEdgesGraphBeforeTruncate = GenerateCollisionEdgeGraph(PLYReader->GetOutput(),collision_edges);
	vtkSmartPointer<vtkMutableUndirectedGraph> collisionEdgesGraphAfterTruncate =  TruncateGraph(collisionEdgesGraphBeforeTruncate);
	
 
	vtkSmartPointer<vtkActor> edge_actor1 = CreateBeforeTruncatePipeline(collisionEdgesGraphBeforeTruncate, collision_edges); 
	vtkSmartPointer<vtkActor> edge_actor2 = CreateAfterTruncatePipeline(collisionEdgesGraphAfterTruncate, collision_edges);
	
	actorEdge1 = edge_actor1;
	actorEdge2 = edge_actor2;
	*/
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(modelReader->GetOutputPort());
 
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetEdgeVisibility(1);
	actor->GetProperty()->SetLineWidth(0.5);


	// Visualize
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
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
	renderer->AddActor(actor);
	renderer->AddActor(actorEdge2);
	
	vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	keypressCallback->SetCallback ( keyPressCallbackFunc );
	renderWindowInteractor->AddObserver ( vtkCommand::KeyPressEvent, keypressCallback );

	vtkSmartPointer<vtkCallbackCommand> pickCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	pickCallback->SetCallback ( pickCallbackFunc );
	pickCallback->SetClientData(TrackballStyle);
	picker->AddObserver(vtkCommand::EndPickEvent,pickCallback);	
	


	ColoredPoint( renderer,modelReader->GetOutput()->GetPoint(sourceVertex), 1.0,0.5,0.0);

	renderWindow->Render(); 	
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
		case '1':
			renderer->RemoveActor(actorEdge1);
			renderer->RemoveActor(actorEdge2);
			renderer->AddActor(actorEdge1);
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			break;
		case '2':
			renderer->RemoveActor(actorEdge1);
			renderer->RemoveActor(actorEdge2);
			renderer->AddActor(actorEdge2);
			renderer->Modified();
			iren->GetRenderWindow()->Render();
			break;
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


bool GenerateGeodesicDistance(int vertexID , std::vector<collisionEdgeInfo>& collision_edges)
{
	if (vertexID < 0 || vertexID >= geosesic_mesh.vertices().size())
		return false;
	geodesic::GeodesicAlgorithmExact exact_algorithm(&geosesic_mesh);
	geodesic::SurfacePoint source(&geosesic_mesh.vertices()[vertexID]);		//create source 
	std::vector<geodesic::SurfacePoint> all_sources(1,source);					//in general, there could be multiple sources, but now we have only one
	exact_algorithm.propagate(all_sources);	//cover the whole mesh

	
	collision_edges.clear();
	for (int i = 0 ; i < exact_algorithm.CollisionEdges().size() ; i++)
	{
		geodesic::edge_pointer e = exact_algorithm.CollisionEdges()[i];	
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
		if (collision_edges[i].type == 1 || collision_edges[i].type == 2)
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
				vtkSmartPointer<vtkAdjacentVertexIterator> adjacent =vtkSmartPointer<vtkAdjacentVertexIterator>::New();
				graph->GetAdjacentVertices(vid[i],adjacent);
				int adjnum = 0;
				while (adjacent->HasNext())
				{
					vtkIdType u = adjacent->Next();
					++adjnum;
				}
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

vtkSmartPointer<vtkActor> CreateBeforeTruncatePipeline(vtkSmartPointer<vtkMutableUndirectedGraph> BTGraph,std::vector<collisionEdgeInfo>& collision_edges)
{
	vtkSmartPointer<vtkGraphToPolyData> g2pBT = vtkSmartPointer<vtkGraphToPolyData>::New(); 
	g2pBT->SetInputData(BTGraph);
	g2pBT->Update();
	vtkSmartPointer<vtkPolyData> collisionEdgesBTPolydata = vtkSmartPointer<vtkPolyData>::New();
	collisionEdgesBTPolydata->ShallowCopy(g2pBT->GetOutput());


	// Setup the colors array
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");
	// Setup two colors 
	unsigned char red[3] = {255, 0, 0};
	unsigned char green[3] = {0, 255, 0};
	unsigned char blue[3] = {0, 0, 255};
	unsigned char yellow[3] = {255, 255, 0};
	int numEdge = collisionEdgesBTPolydata->GetLines()->GetNumberOfCells();
	for (int cellID = 0 ; cellID < collisionEdgesBTPolydata->GetLines()->GetNumberOfCells(); cellID++)
	{
		if (collision_edges[cellID].type == 1)
			colors->InsertNextTupleValue(red);	
		else if (collision_edges[cellID].type == 2)
			colors->InsertNextTupleValue(green);

	}
	
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
	
	collisionEdgesBTPolydata->GetCellData()->SetScalars(colors);

	vtkSmartPointer<vtkPolyDataMapper> edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	edge_mapper->SetInputData(collisionEdgesBTPolydata);
 
	vtkSmartPointer<vtkActor> edge_actor = vtkSmartPointer<vtkActor>::New();
	edge_actor->SetMapper(edge_mapper);
	edge_actor->GetProperty()->SetLineWidth(2.0);

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
	edge_actor->GetProperty()->SetLineWidth(2.0);
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
		pointsActor->GetProperty()->SetPointSize(5);
		pointsActor->GetProperty()->SetColor(r,g,b);

		_renderer->AddActor(pointsActor);
	}
}


