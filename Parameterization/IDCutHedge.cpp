#include "IDCutHedge.h"





IDCutHedge::IDCutHedge():
IDList()
{
	twin = NULL;
	FaceID = -1;
	NewIDSpilt = SPILT_NOT_DETERMINE;
	DuplicateDegree = 0;
	CutDegree = -1;
}

IDCutHedge::IDCutHedge(int dv):
IDList(dv)
{
	twin = NULL;
	FaceID = -1;
	NewIDSpilt = SPILT_NOT_DETERMINE;
	DuplicateDegree = 0;
	CutDegree = -1;
}


IDCutHedge::~IDCutHedge(void)
{
	twin = NULL;
	FaceID = -1;
	NewIDSpilt = SPILT_NOT_DETERMINE;
	CutDegree = -1;
}

void IDCutHedge::AppendVF(int vertexID, int FaceID,int CutDegree ,IDList *tail)
{
	IDCutHedge *now = new IDCutHedge(vertexID);
	now->FaceID = FaceID;
	now->CutDegree = CutDegree;
	IDList *dummy = tail->back;
	now->next = tail;
	tail->back = now;
	now->back = dummy;
	dummy->next =now;
}


 void IDCutHedge::RegisterTwin(IDCutHedge *pCutHedgeH,IDCutHedge *pCutHedgeT)
 {
	if (pCutHedgeH == pCutHedgeT)
		return ;
	else
		pCutHedgeT->ID = pCutHedgeH->next->ID;

	IDCutHedge *now = pCutHedgeH;
	while(nextN(now) != pCutHedgeT)			
	{
		now = (IDCutHedge *)nextN(now);
		if (now->twin == NULL)
		{
			//check for twin
			IDCutHedge *now2 = pCutHedgeT;
			while(backN(now2) != pCutHedgeH)
			{
				now2 = (IDCutHedge *)backN(now2);

				if (now2 == now)
				{
					//not found for sure...
					break;
				}

				if (now2->twin == NULL)
				{
					if (now2->next->ID == now->ID &&
						now2->ID == now->next->ID)
					{
						///found twin!
						now->twin = now2;
						now2->twin = now;
						break;

					}
				}
			}
		}
	}
 }

void IDCutHedge::ReportTwin(IDCutHedge *pCutHedgeH,IDCutHedge *pCutHedgeT,PolarVertex *polarvertex_array)
{
	FILE *fp = NULL;
	fp = fopen("twin_report.txt","wt");
	

	
	IDCutHedge *now = pCutHedgeH;

	while(nextN(now) != pCutHedgeT)			
	{
		now = (IDCutHedge *)nextN(now);
		if (now->twin != NULL)
		{
			//report
			fprintf(fp,"ID:%04d(%f,%f,%f)<--->(%f,%f,%f)%04d:ID\n", now->ID, 
																polarvertex_array[now->ID].p_vertex->x,
																polarvertex_array[now->ID].p_vertex->y,
																polarvertex_array[now->ID].p_vertex->z,

																polarvertex_array[now->twin->ID].p_vertex->x,
																polarvertex_array[now->twin->ID].p_vertex->y,
																polarvertex_array[now->twin->ID].p_vertex->z,
																now->twin->ID);

			
			fprintf(fp,"ID:%04d(%f,%f,%f)<--->(%f,%f,%f)%04d:ID\n", now->next->ID, 
																polarvertex_array[now->next->ID].p_vertex->x,
																polarvertex_array[now->next->ID].p_vertex->y,
																polarvertex_array[now->next->ID].p_vertex->z,

																polarvertex_array[now->twin->next->ID].p_vertex->x,
																polarvertex_array[now->twin->next->ID].p_vertex->y,
																polarvertex_array[now->twin->next->ID].p_vertex->z,
																now->twin->next->ID);
			fprintf(fp,"==========================================\n");
			


		}
	}

	fclose(fp);
}