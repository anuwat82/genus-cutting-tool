#pragma once

#include "IDList.h"
#include "Point3d.h"

class VList: public IDList
{
public:  
  int FaceID; 
  Point3d normal;
  VList();
  VList(int dv);
  VList(int dv,int fid);
  virtual ~VList();
  static void AppendVF(int vertexID, int FaceID, IDList *idList);  

 private:
  VList(const VList& rhs);
  const VList &operator=(const VList& rhs);
};

inline VList * next(VList *p)
{
	return (VList*)p->next;
}
inline VList * back(VList *p)
{
	return (VList*)p->back;
}
