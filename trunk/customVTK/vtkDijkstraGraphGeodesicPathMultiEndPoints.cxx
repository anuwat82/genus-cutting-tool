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
#include "vtkDijkstraGraphGeodesicPathMultiEndPoints.h"

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


vtkStandardNewMacro(vtkDijkstraGraphGeodesicPathMultiEndPoints);

// Construct object with feature angle = 30; all types of edges, except
// manifold edges, are extracted and colored.
vtkDijkstraGraphGeodesicPathMultiEndPoints::vtkDijkstraGraphGeodesicPathMultiEndPoints()
{
	this->EndPointsIdList = NULL;
	EndPointsIdList = vtkIdList::New();
}

vtkDijkstraGraphGeodesicPathMultiEndPoints::~vtkDijkstraGraphGeodesicPathMultiEndPoints()
{
	EndPointsIdList->Delete();
}


//----------------------------------------------------------------------------
int vtkDijkstraGraphGeodesicPathMultiEndPoints::RequestData(
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

  this->ShortestPath( input, this->StartVertex, this->EndPointsIdList );
  this->TraceShortestPath( input, output, this->StartVertex, this->EndVertex );
  return 1;
}


void vtkDijkstraGraphGeodesicPathMultiEndPoints::ShortestPath( vtkDataSet *inData,
                                                int startv, vtkIdList *EndPointsIdList )
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
		if ( u < 0 || u == startv || EndPointsIdList->IsId(u) >= 0 )
          {
          continue;
          }
        this->Internals->BlockedVertices[u] = true;
      }
    }

  this->Internals->CumulativeWeights[startv] = 0;

  this->Internals->HeapInsert(startv);
  this->Internals->OpenVertices[startv] = true;

  bool stop = false;
  while ((u = this->Internals->HeapExtractMin()) >= 0 && !stop)
    {
    // u is now in ClosedVertices since the shortest path to u is determined
    this->Internals->ClosedVertices[u] = true;
    // remove u from OpenVertices
    this->Internals->OpenVertices[u] = false;

    if ( EndPointsIdList->IsId(u) >= 0  && this->StopWhenEndReached)
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

