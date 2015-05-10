#include "VList.h"


VList::VList():
IDList()
{
	FaceID = -1;
	
}

VList::VList(int dv):
IDList(dv)
{
	FaceID = -1;
}

VList::VList(int dv,int fid):
IDList(dv)
{
	FaceID = fid;
}


VList::~VList(void)
{
	/*FaceID = -1;*/
}

void VList::AppendVF(int vertexID, int FaceID,IDList *tail)
{
	VList *now = new VList(vertexID,FaceID);
	now->FaceID = FaceID;	
	IDList *dummy = tail->back;
	now->next = tail;
	tail->back = now;
	now->back = dummy;
	dummy->next =now;
}


