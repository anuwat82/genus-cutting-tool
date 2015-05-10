#pragma once

#include "IDList.h"
#include "MeshStructure.h"

#define SPILT_NOT_DETERMINE	-1	



class IDCutHedge: public IDList
{
public:
  IDCutHedge *twin;
  int FaceID;
  int NewIDSpilt;
  int CutDegree;
  int DuplicateDegree;
  IDCutHedge();
  IDCutHedge(int dv);
  virtual ~IDCutHedge();
  static void AppendVF(int vertexID, int FaceID,int CutDegree, IDList *idList);
  static void RegisterTwin(IDCutHedge *head,IDCutHedge *tail);
  static void ReportTwin(IDCutHedge *head,IDCutHedge *tail, PolarVertex *polarvertex_array);

 private:
  IDCutHedge(const IDCutHedge& rhs);
  const IDCutHedge &operator=(const IDCutHedge& rhs);
};

inline IDCutHedge * next(IDCutHedge *p)
{
	return (IDCutHedge*)p->next;
}
inline IDCutHedge * back(IDCutHedge *p)
{
	return (IDCutHedge*)p->back;
}