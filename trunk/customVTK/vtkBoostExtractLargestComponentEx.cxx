/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkBoostExtractLargestComponent.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkBoostExtractLargestComponentEx.h"

#include "vtkBoostConnectedComponents.h"
#include "vtkDataSetAttributes.h"
#include "vtkDirectedGraph.h"
#include "vtkExecutive.h"
#include "vtkExtractSelectedGraph.h"
#include "vtkGraph.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkSelection.h"
#include "vtkSelectionNode.h"
#include "vtkUndirectedGraph.h"

#include <algorithm>

vtkStandardNewMacro(vtkBoostExtractLargestComponentEx);

vtkBoostExtractLargestComponentEx::vtkBoostExtractLargestComponentEx()
{
	this->OldIdsArray = NULL;
	OldIdsArray = vtkIdTypeArray::New();
}

vtkBoostExtractLargestComponentEx::~vtkBoostExtractLargestComponentEx()
{
	
	OldIdsArray->Delete();
}

int vtkBoostExtractLargestComponentEx::RequestData(vtkInformation *vtkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Get the input and ouptut
  vtkGraph* input = vtkGraph::SafeDownCast(
      inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkGraph* output = vtkGraph::SafeDownCast(
      outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkSmartPointer<vtkGraph> inputCopy;
  if (vtkDirectedGraph::SafeDownCast(input))
    {
    inputCopy = vtkSmartPointer<vtkDirectedGraph>::New();
    }
  else
    {
    inputCopy = vtkSmartPointer<vtkUndirectedGraph>::New();
    }
  inputCopy->ShallowCopy(input);

  // Find all of the connected components
  vtkSmartPointer<vtkBoostConnectedComponents> connectedComponents =
    vtkSmartPointer<vtkBoostConnectedComponents>::New();
  connectedComponents->SetInputData(inputCopy);
  connectedComponents->Update();

  vtkIntArray* components = vtkIntArray::SafeDownCast(
    connectedComponents->GetOutput()->GetVertexData()->GetArray("component"));

  // Create an array to store the count of the number of vertices
  // in every component
  int componentRange[2];
  components->GetValueRange(componentRange);
  std::vector<int> componentCount(componentRange[1] + 1);

  for(vtkIdType i = 0; i < components->GetNumberOfTuples(); i++)
    {
    componentCount[components->GetValue(i)]++;
    }

  // Save the original counts
  std::vector<int> originalComponentCount(componentCount.size());
  std::copy(componentCount.begin(), componentCount.end(), originalComponentCount.begin());

  // Sort in descending order
  std::sort(componentCount.rbegin(), componentCount.rend());

  // Find the original component ID of the component with the highest count
  std::vector<int>::iterator it = find(originalComponentCount.begin(),
    originalComponentCount.end(), componentCount[0]);

  if(it == originalComponentCount.end())
    {
    vtkErrorMacro("Should never get to the end of the components!");
    return 0;
    }

  int largestComponent = it - originalComponentCount.begin();

  vtkDebugMacro("The largest component is " << largestComponent
      << " and it has " << componentCount[0] << " vertices.");


  //store non-iso components number
  NumberOfNonIsoComponents = 0;
  std::vector<int>::iterator it_first1 = find(componentCount.begin(),
    componentCount.end(), 1);
  if (it_first1 != componentCount.end())
  {
	  NumberOfNonIsoComponents = it_first1 - componentCount.begin();
  }

  // Put either the index of the vertices belonging to the largest connected component
  // or the index of the vertices NOT the largest connected component (depending on the
  // InververtSelection flag) into an array to be used to extract this part of the graph.
  vtkSmartPointer<vtkIdTypeArray> ids =
    vtkSmartPointer<vtkIdTypeArray>::New();
  for(vtkIdType i = 0; i < components->GetNumberOfTuples(); i++)
    {
    if(!this->InvertSelection)
      {
      if(components->GetValue(i) == largestComponent)
        {
        ids->InsertNextValue(i);
        }
      }
    else
      {
      if(components->GetValue(i) != largestComponent)
        {
        ids->InsertNextValue(i);
        }
      }
    }

  vtkDebugMacro(<< ids->GetNumberOfTuples() << " values selected.");

  // Mark all of the things in the graph that should be extracted
  vtkSmartPointer<vtkSelection> selection =
      vtkSmartPointer<vtkSelection>::New();

  vtkSmartPointer<vtkSelectionNode> node =
    vtkSmartPointer<vtkSelectionNode>::New();
  selection->AddNode(node);
  node->SetSelectionList(ids);
  node->SetContentType(vtkSelectionNode::INDICES);
  node->SetFieldType(vtkSelectionNode::VERTEX);

  // Extract them
  vtkSmartPointer<vtkExtractSelectedGraph> extractSelectedGraph =
    vtkSmartPointer<vtkExtractSelectedGraph>::New();
  extractSelectedGraph->SetInputData(0, inputCopy);
  extractSelectedGraph->SetInputData(1, selection);
  extractSelectedGraph->Update();

  output->ShallowCopy(extractSelectedGraph->GetOutput());
  
  
  OldIdsArray->DeepCopy(ids);
  return 1;
}

