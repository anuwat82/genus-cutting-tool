/*=========================================================================

  Copyright (c) Karthik Krishnan
  See Copyright.txt for details.

=========================================================================*/

#include "vtkSmartPointer.h"

#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTimerLog.h"
#include "vtkCamera.h"
#include "vtkProperty.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPolyDataNormals.h"
#include "vtkRendererCollection.h"
#include "vtkPolyDataCollection.h"
#include "vtkObjectFactory.h"
#include "vtkIdList.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkContourWidget.h"
#include "vtkOrientedGlyphContourRepresentation.h"
#include "vtkPolygonalSurfacePointPlacer.h"
#include "vtkPolygonalSurfaceContourLineInterpolator2.h"
#include "vtkTestUtilities.h"
#include "vtkTesting.h"

char TestGeodesicPathEventLog[] =
"# StreamVersion 1 i\n"
"RenderEvent 0 0 0 0 0 0 0 i\n"
"EnterEvent 443 322 0 0 0 0 0 i\n"
"LeftButtonPressEvent 304 47 0 0 0 0 0 i\n"
"RenderEvent 304 47 0 0 0 0 0 i\n"
"LeftButtonReleaseEvent 304 47 0 0 0 0 0 i\n"
"RightButtonPressEvent 87 633 0 0 0 0 0 i\n"
"RenderEvent 87 633 0 0 0 0 0 i\n"
"RightButtonReleaseEvent 87 633 0 0 0 0 0 i\n"
"KeyPressEvent 52 503 0 0 113 1 q i\n"
"CharEvent 52 503 0 0 113 1 q i\n"
"ExitEvent 52 503 0 0 113 1 q i\n"
;
int TestGeodesicPath(int argc, char*argv[])
{
  if (argc < 2)
    {
    std::cerr
      << "Demonstrates editing capabilities of a contour widget on polygonal \n"
      << "data.\n"
      << "Usage args: mesh.vtp [Method 0=Dijkstra,1=FastMarching] [InterpolationOrder 0=NearestNeighbor,1=Linear] [height_offset]."
      << std::endl;
    return EXIT_FAILURE;
    }

  vtkNew< vtkXMLPolyDataReader > reader;
  char *fname = vtkTestUtilities::ExpandDataFileName(argc, argv, argv[1]);
  reader->SetFileName(fname);
  reader->Update();
  delete [] fname;

  vtkNew<vtkPolyDataNormals> normals;

  const int geodesicMethod = (argc > 2 ? atoi(argv[2]) : 0);
  const int interpolationOrder = (argc > 3 ? atoi(argv[3]) : 0);
  const double distanceOffset = (argc > 4 ? atof(argv[4]) : 0);

  // We need to ensure that the dataset has normals if a distance offset was
  // specified.
  if (fabs(distanceOffset) > 1e-6)
    {
    normals->SetInputConnection(reader->GetOutputPort());
    normals->SplittingOff();

    // vtkPolygonalSurfacePointPlacer needs cell normals
    // vtkPolygonalSurfaceContourLineInterpolator needs vertex normals
    normals->ComputeCellNormalsOn();
    normals->ComputePointNormalsOn();
    normals->Update();
    }

  vtkPolyData *pd = (fabs(distanceOffset) > 1e-6) ? 
      normals->GetOutput() : reader->GetOutput();

  vtkNew<vtkPolyDataMapper> mapper;
  mapper->SetInputConnection(fabs(distanceOffset) > 1e-6 ? 
      normals->GetOutputPort() : reader->GetOutputPort());

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper.GetPointer());

  vtkNew<vtkRenderer> ren;
  vtkNew<vtkRenderWindow> renWin;
  renWin->AddRenderer(ren.GetPointer());
  vtkNew<vtkRenderWindowInteractor> iren;
  iren->SetRenderWindow(renWin.GetPointer());
   
  ren->AddActor(actor.GetPointer());
  ren->GetActiveCamera()->SetPosition(-3.68, .447, 1.676);
  ren->GetActiveCamera()->Roll(150);
  ren->ResetCamera();
  ren->ResetCameraClippingRange();

  ren->AddActor(actor.GetPointer());

  // Create a contour widget to interactively trace the path

  vtkNew<vtkContourWidget> contourWidget;
  contourWidget->SetInteractor(iren.GetPointer());
  vtkOrientedGlyphContourRepresentation *rep =
    vtkOrientedGlyphContourRepresentation::SafeDownCast(
      contourWidget.GetPointer()->GetRepresentation());
  rep->GetLinesProperty()->SetColor(1, 0.2, 0);
  rep->GetLinesProperty()->SetLineWidth(5.0);

  vtkNew<vtkPolygonalSurfacePointPlacer> pointPlacer;
  pointPlacer->AddProp(actor.GetPointer());
  pointPlacer->GetPolys()->AddItem( pd );
  rep->SetPointPlacer(pointPlacer.GetPointer());

  // Snap the contour nodes to the closest vertices on the mesh
  pointPlacer->SnapToClosestPointOn();

  vtkNew<vtkPolygonalSurfaceContourLineInterpolator2> interpolator;
  interpolator->GetPolys()->AddItem( pd );
  interpolator->SetGeodesicMethod(geodesicMethod);
  interpolator->SetInterpolationOrder(interpolationOrder);
  rep->SetLineInterpolator(interpolator.GetPointer());
  if (fabs(distanceOffset) > 1e-6)
    {
    pointPlacer->SetDistanceOffset( distanceOffset );
    interpolator->SetDistanceOffset( distanceOffset );
    }

  renWin->Render();
  iren->Initialize();
  contourWidget->EnabledOn();

  // Parameters used for the LaTeX figure
  renWin->SetSize(445, 678);
  ren->GetActiveCamera()->SetFocalPoint(0.254678, 0.149946, -0.229227);
  ren->GetActiveCamera()->SetPosition(-0.388447, 0.660936, 0.203816);
  ren->GetActiveCamera()->SetParallelScale(1.06557);
  ren->GetActiveCamera()->SetClippingRange(0.00314092, 3.14092);
  ren->GetActiveCamera()->SetViewUp(-0.182271, -0.759055, 0.62499);

  return vtkTesting::InteractorEventLoop(
      argc, argv, iren.GetPointer(), TestGeodesicPathEventLog );
}
