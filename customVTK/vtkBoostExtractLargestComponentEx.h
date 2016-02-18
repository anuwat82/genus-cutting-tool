/** This is custom vtkBoostExtractLargestComponentEx
   Purpose of this class is to save IDList of input ID.
   Because, vtkBoostExtractLargestComponent will create new IDList and we cannot know how new IDList from vtkBoostExtractLargestComponent represent which IDs in input graph
   And store number of non-isolate component for easier to check.
*/

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
  
  vtkIdTypeArray*  OldIdsArray; //To store ID of input polydata
  int NumberOfNonIsoComponents;
private:
  vtkBoostExtractLargestComponentEx(const vtkBoostExtractLargestComponentEx&);  // Not implemented.
  void operator=(const vtkBoostExtractLargestComponentEx&);  // Not implemented.
};

#endif
