/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkBoostExtractLargestComponent.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkBoostExtractLargestComponent - Extract the largest connected
// component of a graph
//
// .SECTION Description
// vtkBoostExtractLargestComponent finds the largest connected region of a
// vtkGraph. For directed graphs, this returns the largest biconnected component.
// See vtkBoostConnectedComponents for details.

#ifndef vtkBoostExtractLargestComponentEx_h
#define vtkBoostExtractLargestComponentEx_h

#include "vtkInfovisBoostGraphAlgorithmsModule.h" // For export macro
#include "vtkBoostExtractLargestComponent.h"



class VTKINFOVISBOOSTGRAPHALGORITHMS_NO_EXPORT vtkBoostExtractLargestComponentEx : public vtkBoostExtractLargestComponent
{
public:
  vtkTypeMacro(vtkBoostExtractLargestComponentEx, vtkBoostExtractLargestComponent);
  

  // Description:
  // Construct an instance of vtkBoostExtractLargestComponent with
  // InvertSelection set to false.
  static vtkBoostExtractLargestComponentEx* New();

  // Description:
  // Set the flag to determine if the selection should be inverted.
  
  vtkGetMacro(OldIdsArray, vtkIdTypeArray*);
  vtkGetMacro(NumberOfNonIsoComponents, int);
  

protected:
  vtkBoostExtractLargestComponentEx();
  ~vtkBoostExtractLargestComponentEx();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  // Description:
  // Store the choice of whether or not to invert the selection.
  
  vtkIdTypeArray*  OldIdsArray;
  int NumberOfNonIsoComponents;
private:
  vtkBoostExtractLargestComponentEx(const vtkBoostExtractLargestComponentEx&);  // Not implemented.
  void operator=(const vtkBoostExtractLargestComponentEx&);  // Not implemented.
};

#endif
