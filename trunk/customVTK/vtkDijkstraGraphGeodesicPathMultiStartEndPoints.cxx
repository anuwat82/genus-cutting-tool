/*=========================================================================

  Program:   Visualization Toolkit
  Module:  vtkDijkstraGraphGeodesicPath.cxx
  Language:  C++
  Date:    $Date$
  Version:   $Revision$

  Made by Rasmus Paulsen
  email:  rrp(at)imm.dtu.dk
  web:    www.imm.dtu.dk/~rrp/VTK

  This class is not mature enough to enter the official VTK release.
=========================================================================*/
#include "vtkDijkstraGraphGeodesicPathMultiStartEndPoints.h"

#include "vtkCellArray.h"
#include "vtkDijkstraGraphInternals.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"


vtkStandardNewMacro(vtkDijkstraGraphGeodesicPathMultiStartEndPoints);

// Construct object with feature angle = 30; all types of edges, except
// manifold edges, are extracted and colored.
vtkDijkstraGraphGeodesicPathMultiStartEndPoints::vtkDijkstraGraphGeodesicPathMultiStartEndPoints()
{
	//this->EndPointsIdList = NULL;
	//EndPointsIdList = vtkIdList::New();

	//this->StartPointsIdList = NULL;
	//StartPointsIdList = vtkIdList::New();
}

vtkDijkstraGraphGeodesicPathMultiStartEndPoints::~vtkDijkstraGraphGeodesicPathMultiStartEndPoints()
{
	//EndPointsIdList->Delete();
	//StartPointsIdList->Delete();
}


//----------------------------------------------------------------------------
int vtkDijkstraGraphGeodesicPathMultiStartEndPoints::RequestData(
  vtkInformation *           vtkNotUsed( request ),
  vtkInformationVector **    inputVector,
  vtkInformationVector *     outputVector)
{
  vtkInformation * inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo =   outputVector->GetInformationObject(0);

  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  if (!input)
    {
    return 0;
    }

  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  if (!output)
    {
    return 0;
    }

  if ( this->AdjacencyBuildTime.GetMTime() < input->GetMTime() )
    {
    this->Initialize( input );
    }
  else
    {
    this->Reset();
    }

  if (this->NumberOfVertices == 0)
    {
    return 0;
    }

  this->ShortestPath( input, this->StartPointsIdList, this->EndPointsIdList );
  this->TraceShortestPath( input, output, this->StartPointsIdList, this->EndVertex );
  return 1;
}


void vtkDijkstraGraphGeodesicPathMultiStartEndPoints::ShortestPath( vtkDataSet *inData,
                                                 vtkIdList * _startPointsIdList, vtkIdList *_endPointsIdList )
{
  int u, v;

  if( this->RepelPathFromVertices && this->RepelVertices )
    {
    // loop over the pts and if they are in the image
    // get the associated index for that point and mark it as blocked
    for( int i = 0; i < this->RepelVertices->GetNumberOfPoints(); ++i )
      {
        double* pt = this->RepelVertices->GetPoint( i );
        u = inData->FindPoint( pt );
		if ( u < 0 || _startPointsIdList->IsId(u) >=0  || _endPointsIdList->IsId(u) >= 0 )
          {
          continue;
          }
        this->Internals->BlockedVertices[u] = true;
      }
    }
  for (int i = 0 ; i < _startPointsIdList->GetNumberOfIds() ; i++)
  {
	  int  _startV = _startPointsIdList->GetId(i);
	  this->Internals->CumulativeWeights[_startV] = 0;

	  this->Internals->HeapInsert(_startV);
	  this->Internals->OpenVertices[_startV] = true;
  }
  bool stop = false;
  while ((u = this->Internals->HeapExtractMin()) >= 0 && !stop)
    {
    // u is now in ClosedVertices since the shortest path to u is determined
    this->Internals->ClosedVertices[u] = true;
    // remove u from OpenVertices
    this->Internals->OpenVertices[u] = false;

    if ( _endPointsIdList->IsId(u) >= 0  && this->StopWhenEndReached)
      {
      stop = true;
	  this->EndVertex = u; //store endpoint id to EndVertex
      }

    std::map<int,double>::iterator it = this->Internals->Adjacency[u].begin();

    // Update all vertices v adjacent to u
    for ( ; it != this->Internals->Adjacency[u].end(); ++it )
      {
      v = (*it).first;

      // ClosedVertices is the set of vertices with determined shortest path...
      // do not use them again
      if ( !this->Internals->ClosedVertices[v] )
        {
        // Only relax edges where the end is not in ClosedVertices
        // and edge is in OpenVertices
        double w;
        if ( this->Internals->BlockedVertices[v] )
        {
          w = VTK_FLOAT_MAX;
        }
        else
        {
          w = (*it).second + this->CalculateDynamicEdgeCost( inData, u, v );
        }

        if ( this->Internals->OpenVertices[v] )
          {
          this->Relax(u, v, w);
          }
        // add edge v to OpenVertices
        else
          {
          this->Internals->OpenVertices[v] = true;
          this->Internals->CumulativeWeights[v] = this->Internals->CumulativeWeights[u] + w;

          // Set Predecessor of v to be u
          this->Internals->Predecessors[v] = u;
          this->Internals->HeapInsert(v);
          }
        }
      }
    }
}

void vtkDijkstraGraphGeodesicPathMultiStartEndPoints::TraceShortestPath( vtkDataSet* inData, vtkPolyData* outPoly,
                vtkIdList * _startPointsIdList, vtkIdType endv)
{
	{
  vtkPoints   *points = vtkPoints::New();
  vtkCellArray *lines = vtkCellArray::New();

  // n is far to many. Adjusted later
  lines->InsertNextCell(this->NumberOfVertices);

  // trace backward
  int v = endv;
  double pt[3];
  vtkIdType id;
  while ( _startPointsIdList->IsId(v) == -1)
    {
    IdList->InsertNextId(v);

    inData->GetPoint(v,pt);
    id = points->InsertNextPoint(pt);
    lines->InsertCellPoint(id);

    v = this->Internals->Predecessors[v];
    }

  this->IdList->InsertNextId(v);

  inData->GetPoint(v,pt);
  id = points->InsertNextPoint(pt);
  lines->InsertCellPoint(id);

  lines->UpdateCellCount( points->GetNumberOfPoints() );
  outPoly->SetPoints(points);
  points->Delete();
  outPoly->SetLines(lines);
  lines->Delete();
}
}