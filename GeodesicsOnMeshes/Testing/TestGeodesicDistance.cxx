/*=========================================================================

  Copyright (c) Karthik Krishnan
  See Copyright.txt for details.

=========================================================================*/

#include "vtkActor.h"
#include "vtkTestUtilities.h"
#include "vtkRegressionTestImage.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include <assert.h>
#include "vtkLookupTable.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkTimerLog.h"
#include "vtkCamera.h"
#include "vtkPointPicker.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPolyDataNormals.h"
#include "vtkRendererCollection.h"
#include "vtkObjectFactory.h"
#include "vtkFastMarchingGeodesicDistance.h"
#include "vtkIdList.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkScalarBarWidget.h"
#include "vtkScalarBarActor.h"
#include "vtkContourWidget.h"
#include "vtkOrientedGlyphContourRepresentation.h"
#include "vtkPolygonalSurfacePointPlacer.h"
#include "vtkPolygonalSurfaceContourLineInterpolator.h"
#include "vtkProperty.h"
#include "vtkPolyDataCollection.h"

char TestGeodesicDistanceEventLog[] =
"# StreamVersion 1 i\n"
"RenderEvent 0 0 0 0 0 0 0 i\n"
"EnterEvent 446 578 0 0 0 0 0 i\n"
"LeftButtonPressEvent 307 473 0 0 0 0 0 i\n"
"RenderEvent 307 473 0 0 0 0 0 i\n"
"LeftButtonReleaseEvent 307 473 0 0 0 0 0 i\n"
"KeyPressEvent 307 473 0 0 113 1 q i\n"
"CharEvent 307 473 0 0 113 1 q i\n"
"ExitEvent 307 473 0 0 113 1 q i\n"
;

char TestGeodesicDistanceExclusionRegionEventLog[] =
"# StreamVersion 1 i\n"
"RenderEvent 0 0 0 0 0 0 0 i\n"
"EnterEvent 590 131 0 0 0 0 0 i\n"
"LeaveEvent 580 242 0 0 0 0 0 i\n"
"EnterEvent 596 57 0 0 0 0 0 i\n"
"LeftButtonPressEvent 301 322 0 0 0 0 0 i\n"
"RenderEvent 301 322 0 0 0 0 0 i\n"
"LeftButtonReleaseEvent 301 322 0 0 0 0 0 i\n"
"LeftButtonPressEvent 263 276 0 0 0 0 0 i\n"
"RenderEvent 263 276 0 0 0 0 0 i\n"
"LeftButtonReleaseEvent 263 276 0 0 0 0 0 i\n"
"LeftButtonPressEvent 308 225 0 0 0 0 0 i\n"
"RenderEvent 308 225 0 0 0 0 0 i\n"
"LeftButtonReleaseEvent 308 225 0 0 0 0 0 i\n"
"LeftButtonPressEvent 372 275 0 0 0 0 0 i\n"
"RenderEvent 372 275 0 0 0 0 0 i\n"
"LeftButtonReleaseEvent 372 275 0 0 0 0 0 i\n"
"LeftButtonPressEvent 356 322 0 0 0 0 0 i\n"
"RenderEvent 356 322 0 0 0 0 0 i\n"
"LeftButtonReleaseEvent 356 322 0 0 0 0 0 i\n"
"LeftButtonPressEvent 305 322 0 0 0 0 0 i\n"
"RenderEvent 305 322 0 0 0 0 0 i\n"
"RenderEvent 305 322 0 0 0 0 0 i\n"
"LeftButtonReleaseEvent 305 322 0 0 0 0 0 i\n"
"RenderEvent 305 322 0 0 0 0 0 i\n"
"RenderEvent 305 313 0 0 0 0 0 i\n"
"LeftButtonPressEvent 317 295 0 0 0 0 0 i\n"
"RenderEvent 317 295 0 0 0 0 0 i\n"
"LeftButtonReleaseEvent 317 295 0 0 0 0 0 i\n"
"KeyPressEvent 317 295 0 0 113 1 q i\n"
"CharEvent 317 295 0 0 113 1 q i\n"
"ExitEvent 317 295 0 0 113 1 q i\n"
;

char TestGeodesicDistancePropagationWtEventLog[] =
"# StreamVersion 1 i\n"
"RenderEvent 0 0 0 0 0 0 0 i\n"
"EnterEvent 591 27 0 0 0 0 0 i\n"
"LeftButtonPressEvent 328 287 0 0 0 0 0 i\n"
"RenderEvent 328 287 0 0 0 0 0 i\n"
"LeftButtonReleaseEvent 328 287 0 0 0 0 0 i\n"
"KeyPressEvent 328 287 0 0 113 1 q i\n"
"CharEvent 328 287 0 0 113 1 q i\n"
"ExitEvent 328 287 0 0 113 1 q i\n"
;

// This interactor style invokes the fast marching geodesic filter when a
// point is interactively picked.
class FastMarchingHelper : public vtkInteractorStyleTrackballCamera
{
public:
  static FastMarchingHelper* New();
  vtkTypeMacro(FastMarchingHelper, vtkInteractorStyleTrackballCamera);

  virtual void OnLeftButtonDown()
  {
    this->Interactor->GetPicker()->Pick(
      this->Interactor->GetEventPosition()[0],
      this->Interactor->GetEventPosition()[1],
      0, this->Interactor->GetRenderWindow()->
          GetRenderers()->GetFirstRenderer());
    vtkIdType ptId = vtkPointPicker::SafeDownCast(
        this->Interactor->GetPicker())->GetPointId();
    if (ptId >= 0)
      {
      vtkNew< vtkIdList > seeds;
      seeds->InsertNextId(ptId);
      this->Geodesic->SetSeeds(seeds.GetPointer());
      cout << "Seeding using picked pointId: " << ptId << endl;

      // If an exclusion region has been defined, use it.
      if (this->ContourWidget->GetEnabled())
        {
        this->AddExclusionPointIds();
        }

      this->DoFastMarching();
      this->Render();
      }
    else
      {
      vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
      }
  }

  void DoFastMarching()
  {
    cout<<"Computing Geodesic using the fast marching method..."<<endl;
    vtkNew<vtkTimerLog> timer;
    timer->StartTimer();

    this->Geodesic->Update();

    timer->StopTimer();
    cout << "Fast marching took "<<timer->GetElapsedTime() << " s." <<endl;
  }

  void Render()
  {
    this->Mapper->SetInputConnection(this->Geodesic->GetOutputPort());
    this->Mapper->GetLookupTable()->SetRange(0,
        this->Geodesic.GetPointer()->GetMaximumDistance());
    this->Mapper->GetInput()->GetPointData()->SetActiveScalars(
        this->Geodesic.GetPointer()->GetFieldDataName());
    this->Mapper->ScalarVisibilityOn();
    this->ScalarBar->GetScalarBarActor()->SetLookupTable(
      this->Mapper->GetLookupTable());
    this->Mapper->UseLookupTableScalarRangeOn();
    this->ScalarBar->On();
    this->Interactor->Render();
  }

  // Optionally, setup a contour widget, so that the user can interactively
  // trace a closed exclusion region as defined by a contour to confine the
  // fast marching to.
  void SetupExclusionContour(vtkRenderWindowInteractor *iren,
      vtkPolyData *pd, // polydata on which the contour is to be drawn
      vtkActor *surfaceActor /* actor reprsenting this surface */ )
  {
    // Optionally a set of exclusion point ids may be defined to prevent fast
    // marching from bleeding out of certain regions.
    this->ContourWidget->SetRepresentation(this->ContourRep.GetPointer());
    this->ContourWidget->SetInteractor(iren);
    this->ContourRep->GetLinesProperty()->SetColor(1, 1.0, 0);
    this->ContourRep->GetLinesProperty()->SetLineWidth(5.0);

    vtkNew<vtkPolygonalSurfacePointPlacer> pointPlacer;
    pointPlacer->AddProp(surfaceActor);
    pointPlacer->GetPolys()->AddItem( pd );
    pointPlacer->SnapToClosestPointOn();
    this->ContourRep->SetPointPlacer(pointPlacer.GetPointer());

    this->DijkstraInterp->GetPolys()->AddItem( pd );
    this->ContourRep->SetLineInterpolator(this->DijkstraInterp.GetPointer());
    this->ContourWidget->EnabledOn();
  }

  // Optionally add an exclusion region for fast marching
  void AddExclusionPointIds()
  {
    vtkNew< vtkIdList > exclusionPtIds;
    this->DijkstraInterp->GetContourPointIds(
      this->ContourRep.GetPointer(), exclusionPtIds.GetPointer());
    this->Geodesic->SetExclusionPointIds(exclusionPtIds.GetPointer());
  }

  vtkNew< vtkFastMarchingGeodesicDistance > Geodesic;
  vtkNew< vtkPolyDataMapper > Mapper;
  vtkNew< vtkScalarBarWidget > ScalarBar;
  vtkNew< vtkContourWidget > ContourWidget;
  vtkNew< vtkPolygonalSurfaceContourLineInterpolator > DijkstraInterp;
  vtkNew< vtkOrientedGlyphContourRepresentation > ContourRep;
};

vtkStandardNewMacro(FastMarchingHelper);

int TestGeodesicDistance(int argc, char* argv[])
{
  if (argc < 2)
    {
    std::cerr << "Args: Bunny.vtp OutputMesh.vtp [--maxDist distance] "
      << "[--propagationWts PointDataFieldName]" << std::endl;
    }

  vtkNew< vtkXMLPolyDataReader > reader;
  vtkNew<vtkRenderer> ren;
  vtkNew<vtkRenderWindow> renWin;
  vtkNew<vtkActor> actor;
  renWin->AddRenderer(ren.GetPointer());
  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin.GetPointer());

  char *fname = vtkTestUtilities::ExpandDataFileName(argc, argv, argv[1]);
  reader->SetFileName(fname);
  reader->Update();
  delete [] fname;

  vtkNew<FastMarchingHelper> helper;

  char *testStream = TestGeodesicDistanceEventLog;

  // Parse args
  vtkNew< vtkIdList > seeds;
  for (int i = 1; i < argc; ++i)
    {
    if (strcmp(argv[i], "--maxDist") == 0)
      {
      helper->Geodesic->SetDistanceStopCriterion(atof(argv[i+1]));
      }
    if (strcmp(argv[i], "--propagationWts") == 0)
      {
      if (vtkDataArray *propWtArr =
        reader.GetPointer()->GetOutput()->GetPointData()->GetArray(argv[i+1]))
        {
        helper->Geodesic->SetPropagationWeights(propWtArr);
        }
      testStream = TestGeodesicDistancePropagationWtEventLog;
      }
    if (strcmp(argv[i], "--exclusionContour") == 0)
      {
      helper->SetupExclusionContour(iren.GetPointer(),
          reader->GetOutput(), actor.GetPointer());
      testStream = TestGeodesicDistanceExclusionRegionEventLog;
      }
    }

  // Compute normals for shaded surfaces
  vtkNew<vtkPolyDataNormals> normals;
  normals->SetInputConnection(reader.GetPointer()->GetOutputPort());
  normals->SplittingOff();
  normals->ComputeCellNormalsOn();
  normals->ComputePointNormalsOff();
  normals->Update();

  // mapper - actor
  helper->Mapper->SetInputConnection(normals.GetPointer()->GetOutputPort());
  helper->Mapper->ScalarVisibilityOff();

  actor->SetMapper(helper->Mapper.GetPointer());
  ren->AddActor(actor.GetPointer());

  renWin->SetSize(600,600);
  ren->GetActiveCamera()->SetPosition(-3.68, .447, 1.676);
  ren->GetActiveCamera()->Roll(150);
  ren->ResetCamera();
  ren->ResetCameraClippingRange();

  // Create a picker to enable the user to pick a seed interactively
  vtkNew<vtkPointPicker> pointPicker;
  iren->SetPicker(pointPicker.GetPointer());

  // the geodesic filter
  helper->Geodesic->SetInputConnection(reader.GetPointer()->GetOutputPort());

  // Set the point data field name
  helper->Geodesic->SetFieldDataName("FMMDist");

  // Display the distances
  helper->ScalarBar->SetInteractor(iren.GetPointer());
  helper->ScalarBar->GetScalarBarActor()->SetTitle("FMMDist");


  if (seeds->GetNumberOfIds() == 0)
    {
    // interactively pick
    iren->SetInteractorStyle( helper.GetPointer() );
    }
  else // non interactive
    {
    helper->Geodesic->SetSeeds(seeds.GetPointer());
    helper->DoFastMarching();
    helper->Render();
    }

  // Render and interact
  renWin->Render();
  iren->Initialize();

  return vtkTesting::InteractorEventLoop(
    argc, argv, iren.GetPointer(), testStream );
}
