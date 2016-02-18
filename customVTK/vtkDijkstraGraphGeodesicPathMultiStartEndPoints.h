/*=========================================================================
 This is custom vtkDijkstraGraphGeodesicPathMultiStartEndPoints
   Purpose of this class is to support multiple starting vertices and ending vertices.
   vtkDijkstraGraphGeodesicPath (original) supports single source , single target scheme only.
   
*/

#ifndef __vtkDijkstraGraphGeodesicPathMultiStartEndPoints_h
#define __vtkDijkstraGraphGeodesicPathMultiStartEndPoints_h

#include "vtkFiltersModelingModule.h" // For export macro
#include "vtkDijkstraGraphGeodesicPath.h"

class VTKFILTERSMODELING_NO_EXPORT vtkDijkstraGraphGeodesicPathMultiStartEndPoints :
                           public vtkDijkstraGraphGeodesicPath
{
public:
    vtkTypeMacro(vtkDijkstraGraphGeodesicPathMultiStartEndPoints,vtkDijkstraGraphGeodesicPath);
  // Description:
  // Instantiate the class
  static vtkDijkstraGraphGeodesicPathMultiStartEndPoints *New();

  // Description:
  // Standard methods for printing and determining type information.

  
  vtkGetMacro(EndPointsIdList,vtkIdList *);
  vtkSetMacro(EndPointsIdList,vtkIdList *);
  vtkGetMacro(StartPointsIdList,vtkIdList *);
  vtkSetMacro(StartPointsIdList,vtkIdList *);
protected:
  vtkDijkstraGraphGeodesicPathMultiStartEndPoints();
  ~vtkDijkstraGraphGeodesicPathMultiStartEndPoints(); 

  virtual int RequestData(vtkInformation *, vtkInformationVector **,
                          vtkInformationVector *);

  // Calculate shortest path from vertex startv to vertex endv.
  virtual void ShortestPath( vtkDataSet *inData, vtkIdList * _startPointsIdList, vtkIdList *_endPointsIdList );

  virtual void TraceShortestPath( vtkDataSet* inData, vtkPolyData* outPoly,
                vtkIdList * _startPointsIdList, vtkIdType endv);
  
  vtkIdList *EndPointsIdList;  //target vertices ID  based on input polygondata
  vtkIdList *StartPointsIdList; //source vertices ID based on input polygondata
private:
  vtkDijkstraGraphGeodesicPathMultiStartEndPoints(const vtkDijkstraGraphGeodesicPathMultiStartEndPoints&);  // Not implemented.
  void operator=(const vtkDijkstraGraphGeodesicPathMultiStartEndPoints&);  // Not implemented.

};

#endif

