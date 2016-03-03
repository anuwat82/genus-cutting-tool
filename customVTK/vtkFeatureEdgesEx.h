/** This is custom vtkFeatureEdgeEx
   Purpose of this class is to save IDList of input ID.
   Because, vtkFeatureEdges will create new IDList and we cannot know how new IDList from vtkFeatureEdges represent which IDs in input PolyData
*/


#ifndef __vtkFeatureEdgesEx_h
#define __vtkFeatureEdgesEx_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkFeatureEdges.h"

class VTKFILTERSCORE_EXPORT vtkFeatureEdgesEx : public vtkFeatureEdges
{
public:
  vtkTypeMacro(vtkFeatureEdgesEx,vtkFeatureEdges);
  static vtkFeatureEdgesEx *New();
  vtkIdType GetOldIdFromCurrentID(vtkIdType currentID);
  vtkIdList * GetOldIdList() {return oldIdList; }
protected:
  vtkFeatureEdgesEx();
  ~vtkFeatureEdgesEx();
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  

  vtkIdList *oldIdList; //To store ID of input polydata
private:
  vtkFeatureEdgesEx(const vtkFeatureEdgesEx&);  // Not implemented.
  void operator=(const vtkFeatureEdgesEx&);  // Not implemented.
};

#endif


