#include "stdafx.h"
#include <Windows.h>
#include "utils.h"
#include "VTK_Header.h"
#include "geodesicVTK.h"
#include "Parameterization\PolygonsData.h"

void keyPressCallbackFunc(vtkObject*, unsigned long eid, void* clientdata, void *calldata);
void pickCallbackFunc(vtkObject*, unsigned long eid, void* clientdata, void *calldata);

vtkWeakPointer<vtkTexture> checkboard_texture;
vtkWeakPointer<vtkTexture> image_texture;
int main(int argc, char* argv[])
{
	vtkObject::GlobalWarningDisplayOff();

	vtkSmartPointer<vtkPNGReader> pngReader = vtkSmartPointer<vtkPNGReader>::New();
	pngReader->SetFileName("./texture/square_texture.png");
	pngReader->Update();

	vtkSmartPointer<vtkJPEGReader> jpgReader = vtkSmartPointer<vtkJPEGReader>::New();
	jpgReader->SetFileName("./texture/check_board.jpg");
	jpgReader->Update();
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
		if (!GetFileName(filename,__T("Open model file..."),"PLY Files\0*.ply\0OFF files\0*.off\0All Files\0*.*\0"))
		{
			std::cout << "Exit due to no model file input." << endl;
			return -1;
		}
	}
	std::string ext = GetFileExtension(filename);
	vtkSmartPointer<vtkPolyDataAlgorithm> modelReader;

	
	vtkSmartPointer<vtkTexture> square_texture = vtkSmartPointer<vtkTexture>::New();
	square_texture->SetInputConnection(pngReader->GetOutputPort());
	square_texture->SetBlendingMode(vtkTexture::VTKTextureBlendingMode::VTK_TEXTURE_BLENDING_MODE_ADD  );
	square_texture->PremultipliedAlphaOn(); 	
	image_texture = square_texture;
	vtkSmartPointer<vtkTexture> _checkboard_texture = vtkSmartPointer<vtkTexture>::New();
	_checkboard_texture->SetInputConnection(jpgReader->GetOutputPort());
	_checkboard_texture->SetBlendingMode(vtkTexture::VTKTextureBlendingMode::VTK_TEXTURE_BLENDING_MODE_REPLACE  );
	_checkboard_texture->PremultipliedAlphaOn(); 	
	checkboard_texture = _checkboard_texture;
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
	{
		cout << "************************************" << endl;
		cout << "Loaded file: " << filename << endl;
		cout << "Number Of Vertices: " << polydata->GetNumberOfPoints() << "\n";
		cout << "Number Of Lines: " << polydata->GetNumberOfLines() << "\n";
		cout << "Number Of Polygons: " << polydata->GetNumberOfPolys() << "\n";
		cout << "Number Of Triangle Strips: " << polydata->GetNumberOfStrips() << "\n";

		cout << "Number Of Pieces: " << polydata->GetNumberOfPieces() << endl;
		cout << "Piece: " << polydata->GetPiece() << endl;
		cout << "Ghost Level: " << polydata->GetGhostLevel() << endl;
		cout << "************************************" << endl;
	}

	//Process(modelReader->GetOutput() , sourceVertex);


	// Create the tree
	vtkSmartPointer<vtkOBBTree> obbTree = vtkSmartPointer<vtkOBBTree>::New();  
	obbTree->SetDataSet(polydata);
	obbTree->BuildLocator();
	
	

	
	vtkSmartPointer<vtkPolyDataMapper> _mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	_mapper->SetInputConnection(modelReader->GetOutputPort());
	mapper = _mapper;
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
	
	actorMainPoly->SetTexture(checkboard_texture);
	
	
	vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	keypressCallback->SetCallback ( keyPressCallbackFunc );
	renderWindowInteractor->AddObserver ( vtkCommand::KeyPressEvent, keypressCallback );
	

	vtkSmartPointer<vtkCallbackCommand> pickCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	pickCallback->SetCallback ( pickCallbackFunc );
	pickCallback->SetClientData(TrackballStyle);
	picker->AddObserver(vtkCommand::EndPickEvent,pickCallback);	
	/*
	{
		vtkCamera *cam  = renderer->GetActiveCamera();
		double pos[3] = {7.21609, -25.5471, 20.343};
		
		double up[3] = {-0.17565, 0.582402, 0.793697};
		
		cam->SetPosition(pos);
		cam->SetViewUp(up);
		
	}
	*/
	 ColoredPoint( renderer,modelReader->GetOutput()->GetPoint(sourceVertex), 1.0,0.5,0.0);
	 /*
	 vtkSmartPointer<vtkFeatureEdges> featureEdges =
    vtkSmartPointer<vtkFeatureEdges>::New();
	 featureEdges->SetInputData(polydata);
  featureEdges->BoundaryEdgesOn();
  featureEdges->FeatureEdgesOff();
  featureEdges->ManifoldEdgesOff();
  featureEdges->NonManifoldEdgesOff();
  featureEdges->Update();
  vtkSmartPointer<vtkPolyDataMapper> bedgeMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  bedgeMapper->SetInputConnection(featureEdges->GetOutputPort());
  vtkSmartPointer<vtkActor> bedgeActor =
    vtkSmartPointer<vtkActor>::New();
  bedgeActor->SetMapper(bedgeMapper);
  bedgeActor->GetProperty()->SetLineWidth(2);
	//ColorMesh();
	*/
	renderer->AddActor(actorEdge2);
	renderer->AddActor(actorPoly1);
	//renderer->AddActor(bedgeActor);
	
	renderer->SetBackground(.2, .2, .2);
	//renderer->ResetCamera(-10,10,-10,10,-10,10);
	renderWindow->Render(); 	

	cout << "===============================" << endl <<
		"CTRL + Left click to change starting point" <<endl <<
		"Function Key:" << endl <<
		"F6 = homotopy cutting" << endl <<
			"  1 = display both kind of crossing edges" << endl <<
			"  2 = display edges after removing dangling ones" << endl <<
			"  3 = display our final cutting path (red)" << endl <<
			"  4 = display original final cutting path (yellow)" << endl <<
		"F7 = iterative augmented cutting" << endl <<
		"F8 = square parameterizations (25% brute force)" << endl <<
		"F9 = square parameterizations (step sampling)" << endl << endl <<
		
		"F10 = save screenshot" << endl << 
		"+/- = increase/decrease opacity" << endl << 

		"===============================" << endl;


	renderWindowInteractor->Start();
	

	if (exact_algorithm)
	{
		delete exact_algorithm;
		exact_algorithm = NULL;
	}


	return 0;
}



void pickCallbackFunc(vtkObject* caller, unsigned long eid, void* clientdata, void *calldata)
{
	
	vtkPointPicker* picker = static_cast<vtkPointPicker*>(caller);
	MouseInteractorStylePP* TrackballStyle = static_cast<MouseInteractorStylePP*>(clientdata);

	if (picker->GetPointId() >= 0 )
	{
		int vertexID = picker->GetPointId();
		

		if (TrackballStyle->GetInteractor()->GetControlKey() != 0)
		{
			vtkSmartPointer<vtkPolyData> PolyData = polydata;
			InitialGeodesic(PolyData,vertexID);
			ColoredPoint( TrackballStyle->GetInteractor()->GetRenderWindow()->GetRenderers()->GetFirstRenderer(),PolyData->GetPoint(vertexID), 0.0,1.0,0.0);
			actorEdge1->Modified();
			actorEdge2->Modified();
			TrackballStyle->GetInteractor()->GetRenderWindow()->Render();
			std::cout << "Picked PointID:" << vertexID <<  std::endl;
		}
	}
}

void keyPressCallbackFunc(vtkObject* caller, unsigned long eid, void* clientdata, void *calldata)
{
	vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
	vtkRenderer *renderer = static_cast<vtkRenderer*>(iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer());
 
	
	char *ch = iren->GetKeySym();
	if (strlen(ch) == 1)
	{
	switch (*ch)
	{
		case 'r':
			{
				renderer->ResetCamera();
				renderer->Modified();
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
				if (exact_algorithm == NULL)	
					InitialGeodesic(polydata.GetPointer(),sourceVertex);
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
			/*
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
			*/
		case 'm':
			{
			vtkTexture *currentTexture = actorPoly1->GetTexture();
			if (currentTexture == NULL)
				actorPoly1->SetTexture(checkboard_texture);
			else if (currentTexture == checkboard_texture)
				actorPoly1->SetTexture(image_texture);
			else if (currentTexture == image_texture)
				actorPoly1->SetTexture(NULL);
			}
			break;
		case 'c':
			//cutting 
			if (modTruncate.isReadyToCut())
			{
				std::string filename;
				if (GetFileName(filename,__T("Save disk topology model file..."),"PLY File\0*.ply\0All Files\0*.*\0",true))
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
	}	
	std::string key(ch);
	if (key == "F10")
	{
		ScreenShot(iren->GetRenderWindow());
	}
	else if (key == "F6")
	{
		//homotopy cutting		
		if (exact_algorithm == NULL)	
			InitialGeodesic(polydata.GetPointer(),sourceVertex);
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
		cout.precision(15);
		cout << "Homotopy Cutting finished..."<< endl;
		printf("time consume original: %.07f sec\n", time_original);// << "sec" << endl;
		printf("time consume proposed: %.07f sec\n", time_proposed);// << "sec" << endl;
		//cout << "time consume proposed:" << time_proposed << "sec" << endl;
		cout << "========================================" << endl;
		printf("Do not forget to call iterated augment cutting (F7) before perform any parameterization!\n");
	}
	else if (key == "F7")
	{
		//iterated augment cutting
		
		vtkSmartPointer<vtkPolyData> outputOptimumOriginal = vtkSmartPointer<vtkPolyData>::New();
		double time_original = 0;
		
		if (originalTruncate.isReadyToCut())
		{
			cout << "Original Iterated Augment Cutting Started..."<< endl;
			vtkSmartPointer<vtkPolyData> diskOriginal = originalTruncate.GetDiskTopologyPolydata();
			CPolygonsData polygon ;
			polygon.InitailDiskTopology(diskOriginal);
			polygon.IteratedAugmentCutOriginal(&time_original,outputOptimumOriginal.GetPointer());
			cout << "Original Iterated Augment Cutting Finished..."<< endl;

			//adjust edge
			vtkSmartPointer<vtkMutableUndirectedGraph> newGraph = CreateBoundaryGraph(outputOptimumOriginal);
			
			vtkSmartPointer<vtkGraphToPolyData> g2pAT = vtkSmartPointer<vtkGraphToPolyData>::New(); 
			g2pAT->SetInputData(newGraph);
			g2pAT->Update();

			vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
			colors->SetNumberOfComponents(3);
			colors->SetName("Colors");
			// Setup two colors 
			unsigned char red[3] = {255, 0, 0};
			unsigned char yellow[3] = {255, 255, 0};
			int numEdge = g2pAT->GetOutput()->GetLines()->GetNumberOfCells();
			for (int cellID = 0 ; cellID < numEdge; cellID++)
			{
				colors->InsertNextTupleValue(yellow);		
			}
			g2pAT->GetOutput()->GetCellData()->SetScalars(colors);

			vtkSmartPointer<vtkPolyDataMapper> __edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			__edge_mapper->SetInputData(g2pAT->GetOutput());
			actorEdge4->SetMapper(__edge_mapper);
			


		}
		
		double time_proposed = 0;
		vtkSmartPointer<vtkPolyData> outputOptimumPropose = vtkSmartPointer<vtkPolyData>::New();
		if (modTruncate.isReadyToCut())
		{
			cout << "Proposed Iterated Augment Cutting Started..."<< endl;
			vtkSmartPointer<vtkPolyData> diskPropose = modTruncate.GetDiskTopologyPolydata();
			CPolygonsData polygon ;
			
			polygon.InitailDiskTopology(diskPropose);
			polygon.IteratedAugmentCut(&time_proposed,outputOptimumPropose.GetPointer());

			cout << "Proposed Iterated Augment Cutting Finished..."<< endl;

			vtkSmartPointer<vtkMutableUndirectedGraph> newGraph = CreateBoundaryGraph(outputOptimumPropose);
			vtkSmartPointer<vtkGraphToPolyData> g2pAT = vtkSmartPointer<vtkGraphToPolyData>::New(); 
			g2pAT->SetInputData(newGraph);
			g2pAT->Update();

			vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
			colors->SetNumberOfComponents(3);
			colors->SetName("Colors");
			// Setup two colors 
			unsigned char red[3] = {255, 0, 0};
			unsigned char yellow[3] = {255, 255, 0};
			int numEdge = g2pAT->GetOutput()->GetLines()->GetNumberOfCells();
			for (int cellID = 0 ; cellID < numEdge; cellID++)
			{
				colors->InsertNextTupleValue(red);		
			}
			g2pAT->GetOutput()->GetCellData()->SetScalars(colors);


			vtkSmartPointer<vtkPolyDataMapper> __edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			__edge_mapper->SetInputData(g2pAT->GetOutput());
			actorEdge3->SetMapper(__edge_mapper);


			disk_polydata = vtkSmartPointer<vtkPolyData>::New();
			disk_polydata->ShallowCopy(outputOptimumPropose);
		}

		

		cout << "Iterated Augment Cuttings finished..."<< endl;
		cout << "time consume original:" << time_original << "sec" << endl;
		cout << "time consume proposed:" << time_proposed << "sec" << endl;
		cout << "========================================" << endl;

		char answer;
		cout << "Do you want to save optimum disk topology meshes? (y/n):";
		cin >> answer;

		if (answer == 'y' ||answer == 'Y')
		{
			std::string filename;
			if (GetFileName(filename,__T("Save optimal disk topology model file..."),"PLY File\0*.ply\0All Files\0*.*\0",true))
			{
				
				vtkSmartPointer<vtkPLYWriter> plyWriter = vtkSmartPointer<vtkPLYWriter>::New();
					
				std::string propose_filename = filename + std::string("_propose.ply");
					std::string original_filename = filename +  std::string("_original.ply");
					
				plyWriter->SetInputData(outputOptimumPropose);
				plyWriter->SetFileName(propose_filename.c_str());
				plyWriter->Update();
				cout << "Disk topology mesh export :" << propose_filename.c_str() << endl;

				plyWriter->SetInputData(outputOptimumOriginal);
				plyWriter->SetFileName(original_filename.c_str());
				plyWriter->Update();					
				cout << "Disk topology mesh export :" << original_filename.c_str() << endl;
				
			}
		}
		else
		{
			cout << "No save" << endl;
		}

		
	}
	else if (key == "F8")
	{
		//Sqaure Parameterization brute force

		vtkSmartPointer<vtkPolyData> inputPolydata;

		vtkSmartPointer<vtkFloatArray> texCoord = vtkSmartPointer<vtkFloatArray>::New();
		if (disk_polydata.GetPointer() != NULL && disk_polydata->GetNumberOfPoints() > 0)
		{
			//use  disk_polydata as input
			vtkSmartPointer<vtkPolyDataMapper> diskPolyDataMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			diskPolyDataMapper->SetInputData(disk_polydata);
			actorMainPoly->SetMapper(diskPolyDataMapper);
			inputPolydata = disk_polydata;		
			

		}
		else
		{
			//use polydata as input (load disk topolopy file directly)
			inputPolydata = polydata;			
		}
		double calTime;
		unsigned int calCount;
		char answer;
		cout << "Do you want to manually set test case? (y/n):";
		cin >> answer;

		if (answer == 'y' ||answer == 'Y')
		{

			int manual_case;
			do
			{			
				cout << endl << "Enter test case value:";
				while(!(cin >> manual_case)){
					cout << "Bad value!" << endl;
					cout << "Re enter the value:";
					cin.clear();
					cin.ignore(numeric_limits<streamsize>::max(), '\n');
				}
			}
			while (manual_case < 0 );
			CPolygonsData polygon ;
			polygon.InitailDiskTopology(inputPolydata);			
			polygon.SquareParameterizationManual(manual_case,&calTime,texCoord.GetPointer());
		}
		else
		{
			cout << "Perform Parameterization 25% brute force..."<< endl;
			CPolygonsData polygon ;
			polygon.InitailDiskTopology(inputPolydata);
			
		
			polygon.SquareParameterizationOptimization(1,&calCount,&calTime,texCoord.GetPointer());
		}
		inputPolydata->GetPointData()->SetTCoords(texCoord);
		cout << "time consume: " << calTime << " sec" << endl;
		cout << "total test cases examined: " << calCount << " times" << endl;
		cout << "========================================" << endl;
		AskForSaveSqp(texCoord);
		//AskForSaveParameterizationPLY(inputPolydata,texCoord);
		
	}
	else if (key == "F9")
	{
		//Sqaure Parameterization step sampling
		
		vtkSmartPointer<vtkPolyData> inputPolydata;
		if (disk_polydata.GetPointer() != NULL && disk_polydata->GetNumberOfPoints() > 0)
		{
			//use  disk_polydata as input
			vtkSmartPointer<vtkPolyDataMapper> diskPolyDataMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			diskPolyDataMapper->SetInputData(disk_polydata);
			actorMainPoly->SetMapper(diskPolyDataMapper);
			inputPolydata = disk_polydata;	
		}
		else
		{
			//use polydata as input (load disk topolopy file directly)
			inputPolydata = polydata;			
		}

		unsigned int step_value = 0;  // 0 means using formula
		char answer;
		cout << "Do you want to manually set step value? (y/n):";
		cin >> answer;

		if (answer == 'y' ||answer == 'Y')
		{
			vtkSmartPointer<vtkFeatureEdges> borderEdges = vtkSmartPointer<vtkFeatureEdges>::New();
			borderEdges->SetInputData(inputPolydata);
			borderEdges->SetBoundaryEdges(1);
			borderEdges->SetFeatureEdges(0);
			borderEdges->SetNonManifoldEdges(0);
			borderEdges->Update();
			int numBedge = borderEdges->GetOutput()->GetNumberOfLines();
			unsigned int limit = (unsigned int )(numBedge*0.25*0.5);
			do
			{			
				cout << endl << "Enter step value (2 ~ " << limit << "):";
				while(!(cin >> step_value)){
					cout << "Bad value!" << endl;
					cout << "Re enter the value:";
					cin.clear();
					cin.ignore(numeric_limits<streamsize>::max(), '\n');
				}
			}
			while (step_value < 2 || step_value > limit);
			
			cout << "Use step value " << step_value << endl;
		}
		else
		{

			cout << "Use formula step value" << endl;
		}


		CPolygonsData polygon ;
		polygon.InitailDiskTopology(inputPolydata);
		double calTime;
		unsigned int calCount;
		vtkSmartPointer<vtkFloatArray> texCoord = vtkSmartPointer<vtkFloatArray>::New();
		polygon.SquareParameterizationOptimization(step_value,&calCount,&calTime,texCoord.GetPointer());
		inputPolydata->GetPointData()->SetTCoords(texCoord);
		
		cout << "consuming time : " << calTime << " sec" << endl;
		cout << "total test cases examined: " << calCount << " times" << endl;
		cout << "========================================" << endl;
		AskForSaveSqp(texCoord);
		AskForSaveParameterizationPLY(inputPolydata,texCoord);	
	}
	else if (key == "F3")
	{
		
		cout << "Perform Circular Parameterization ..."<< endl;
		vtkSmartPointer<vtkPolyData> inputPolydata;
		if (disk_polydata.GetPointer() != NULL && disk_polydata->GetNumberOfPoints() > 0)
		{
			//use  disk_polydata as input
			vtkSmartPointer<vtkPolyDataMapper> diskPolyDataMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			diskPolyDataMapper->SetInputData(disk_polydata);
			actorMainPoly->SetMapper(diskPolyDataMapper);
			inputPolydata = disk_polydata;		
		}
		else
		{
			//use polydata as input (load disk topolopy file directly)
			inputPolydata = polydata;			
		}
		CPolygonsData polygon ;
		polygon.InitailDiskTopology(inputPolydata);
		double calTime;
		unsigned int calCount;
		vtkSmartPointer<vtkDoubleArray> stretch = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkFloatArray> texCoord = vtkSmartPointer<vtkFloatArray>::New();
		//polygon.CircleParameterizationOptimization(&calTime,NULL,stretch.GetPointer());
		//ColorMeshFace(stretch);
		//polygon.CheckBoundaryMapping(NULL);
		polygon.SquareParameterizationExperiment(&calTime,texCoord.GetPointer());
		inputPolydata->GetPointData()->SetTCoords(texCoord);
		//ColorMeshFace(stretch);
		cout << "time consume: " << calTime << " sec" << endl;		
		
		cout << "========================================" << endl;
		AskForSaveParameterizationPLY(inputPolydata,texCoord);
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

