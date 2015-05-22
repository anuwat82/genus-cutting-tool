#define BOOST_LIB_DIAGNOSTIC
#include "MyParameterization.h"
#include "GPUSolver.cuh"
//#include <NL/nl.h>

//#include <boost/thread.hpp>
#include <ppl.h>
bool IsNumber(double x);

MyParameterization::MyParameterization(void)
{

	// Default setting
    pickID=-1;
    paramtype=2;
    boundarytype=0;
	//boundarytype=0;
    weighttype=0;
    iteNum=2000;
    gammaP=0.75;
    PCBCGerror=pow(0.1,6.0);
    smooth=1;
    intrinsiclambda=0.5;
    boundarysigma=1;
	mirrorFace = NULL;

	gammaU = 1.0;
	gammaV = 1.0;
	
	PT = NULL;
	sigmaU = NULL;
	sigmaV = NULL;
	allocateVsize = allocateFsize = 0;

	IDtool = new IDSet();
	PT = new PointTool();
}


MyParameterization::~MyParameterization(void)
{
	if (mirrorFace!= NULL)
	{
		delete [] mirrorFace;
		mirrorFace = NULL;
	}
	if (sigmaU != NULL)
		delete [] sigmaU;
	if (sigmaV != NULL)
		delete [] sigmaV;
}

void	MyParameterization::DeleteFace(int id)
{
	if (Face[id])
	{
		int di = Face[id][0];
		int dj = Face[id][1];
		int dk = Face[id][2];

		IDtool->DeleteF(id,FHead[di],FTail[di]);
		IDtool->DeleteF(id,FHead[dj],FTail[dj]);
		IDtool->DeleteF(id,FHead[dk],FTail[dk]);

  
		/* One */
		neighborF[di]--;
  
		//must delete VF before I because to check edge that have half-edge or not.
		IDtool->DeleteVF(dj,dk,VHead[di],VTail[di]);
		IDtool->DeleteI(dj,IHead[di],ITail[di],VHead[di],VTail[di],di,neighborI);
		IDtool->DeleteI(dk,IHead[di],ITail[di],VHead[di],VTail[di],di,neighborI);

		/* Two */
		neighborF[dj]--;
  
		IDtool->DeleteVF(dk,di,VHead[dj],VTail[dj]);
		IDtool->DeleteI(di,IHead[dj],ITail[dj],VHead[dj],VTail[dj],dj,neighborI);
		IDtool->DeleteI(dk,IHead[dj],ITail[dj],VHead[dj],VTail[dj],dj,neighborI);		
		
		/* Three */
		neighborF[dk]--;
  
		IDtool->DeleteVF(di,dj,VHead[dk],VTail[dk]);
		IDtool->DeleteI(di,IHead[dk],ITail[dk],VHead[dk],VTail[dk],dk,neighborI);		
		IDtool->DeleteI(dj,IHead[dk],ITail[dk],VHead[dk],VTail[dk],dk,neighborI);		
		
		delete Face[id];
		Face[id] = NULL;
	}
}

void    MyParameterization::SetPolarVertexAndFaceAndBorder(	PolarVertex *pIPV,
															int *pIOnum_PV,
															int num_originalVertex,
															int	*pIOFace,
															int *pIOnum_Face,
															IDCutHedge *pCutHedgeH,
															IDCutHedge *pCutHedgeT,
															int num_valen2Boundary,
															bool handleValence2)
{
	//to prevent degeneracy we cannot allow boundary vertex to be valence-2.
	
	//max possible  num face is (*pIOnum_Face)+4
	//max possible  num polar vertex is num_PV+2
	if (num_valen2Boundary > 0)
	{
		memoryallocate( (*pIOnum_PV)+num_valen2Boundary, (*pIOnum_Face)+(num_valen2Boundary*2));
	}
	else
	{
		memoryallocate( (*pIOnum_PV)+2, (*pIOnum_Face)+4);
	}
	//sigmaU = new double[numberV];
	//sigmaV = new double[numberV];
	numberV = (*pIOnum_PV);
	numberF = (*pIOnum_Face);

#pragma omp parallel for
	for(int i=0;i<num_originalVertex;i++)
	{		
		setPoint(i,pIPV[i].p_vertex->x,pIPV[i].p_vertex->y,pIPV[i].p_vertex->z);
	}
	
//#pragma omp parallel for	
	for(int i=0;i<numberF;i++)
	{		
		setFace(i,pIOFace[(i*3)+0],pIOFace[(i*3)+1],pIOFace[(i*3)+2]);
	}
	

	//to do: how to store twin halfedge
	//re-structure cut-route.	
	//boundarySurfaceFaceInfo = 
	if (pCutHedgeH)
	{
		int vertexNewIndex = num_originalVertex;
		IDCutHedge *now = pCutHedgeH;
		while(next(now) != pCutHedgeT)			
		{
		
			now = next(now);		
			if (now->DuplicateDegree > 0)
			{
				//create new point
				pIPV[vertexNewIndex].p_vertex = pIPV[now->ID].p_vertex;
				now->NewIDSpilt = vertexNewIndex;
				setPoint(vertexNewIndex,pIPV[vertexNewIndex].p_vertex->x,pIPV[vertexNewIndex].p_vertex->y,pIPV[vertexNewIndex].p_vertex->z);
				vertexNewIndex++;

				//find starthalfege;
				int startHalfEdge = -1;
				for (int i = 0 ; i < 3 ; i++)
				{
					if (Face[now->FaceID][i] != now->ID && Face[now->FaceID][i] != now->next->ID)
					{
						startHalfEdge = Face[now->FaceID][i];
						break;
					}
				}

				//determine stophalfEdge
				int stopHalfEdge = back(now)->NewIDSpilt;		
				ChangeVertexIndexInSurfaceFaces(now->ID, now->NewIDSpilt, now->FaceID,startHalfEdge , stopHalfEdge,pIOFace);
				for (int faceIdx = 0 ; faceIdx < 3 ; faceIdx++)
				{
					if (pIOFace[(now->FaceID*3) + faceIdx]  == now->ID)
					{
						pIOFace[(now->FaceID*3) + faceIdx] = now->NewIDSpilt;
						break;
					}				
				}
				DeleteFace(now->FaceID);
				setFace(now->FaceID,pIOFace[(now->FaceID*3)+0],pIOFace[(now->FaceID*3)+1],pIOFace[(now->FaceID*3)+2]);
			
			}
			else
			{
				now->NewIDSpilt = now->ID; //use same id  no change
			}
		}
	}
	SetBoundaryLines();
	
	HandleValence2Boundary(	
							pIPV,
							pIOnum_PV,
							pIOFace,
							pIOnum_Face,
							pCutHedgeH,
							pCutHedgeT
							);
	
	IDCutHedge::RegisterTwin(pCutHedgeH,pCutHedgeT);


	SortV();
	setAreaMap3D();
	CalculateEdgeLength();
	//createMirrorFaceData();
	
}
bool	MyParameterization::HandleValence2Boundary(
													PolarVertex *pIPV,
													int *pIOnum_PV,
													int	*pIOFace,
													int *pIOnum_Face,													
													IDCutHedge *pCutHedgeH,
													IDCutHedge *pCutHedgeT)
{

	bool change = false;
	for(int i=0;i<numberV;i++)
	{
		if(boundary[i]==1 && neighborF[i] == 1)
		{
			//divide that face into 2 face.
			change = true;
			
			//create new vertex , polar ... center of edge
			int v1 = IHead[i]->next->ID;
			int v2 = ITail[i]->back->ID;
			int f0 = FHead[i]->next->ID;		
			Vertex *newPT = new Vertex;
			
			newPT->x = (point[v1]->x + point[v2]->x)/2.0;
			newPT->y = (point[v1]->y + point[v2]->y)/2.0;
			newPT->z = (point[v1]->z + point[v2]->z)/2.0;
		

			//newCreateVertex.push_back(newPT);
			pIPV[(*pIOnum_PV)].p_vertex =newPT;// &(newCreateVertex.back());
			int newPT_ID = (*pIOnum_PV);
			setPoint(newPT_ID,newPT->x,newPT->y,newPT->z);
			(*pIOnum_PV)++;
			numberV++;
			
			
			//find another face that share edge v1-v2 with face f0
			int f1 = -1;
			int i1 = -1;
			IDList *nowf = FHead[v1];
			bool found = false;
			while (nextN(nowf) != FTail[v1] && !found)
			{
				nowf = nextN(nowf);
				if (nowf->ID != f0)
				{
					for (int j = 0 ; j < 3 ; j++)
					{
						if (Face[nowf->ID][j] == v2)
						{
							f1 = nowf->ID;
							if (Face[nowf->ID][(j+1)%3] != v1)
							{
								i1 = Face[nowf->ID][(j+1)%3];
							}
							else
							{
								i1 = Face[nowf->ID][(j+2)%3];
							}
							found = true;
							break;
						}
					}
				}
			}
			DeleteFace(f0);
			DeleteFace(f1);
		
			setFace(f0,newPT_ID,i,v1);
			setFace(f1,newPT_ID,v1,i1);
			pIOFace[f0*3 + 0] = newPT_ID;
			pIOFace[f0*3 + 1] = i;
			pIOFace[f0*3 + 2] = v1;
			pIOFace[f1*3 + 0] = newPT_ID;
			pIOFace[f1*3 + 1] = v1;
			pIOFace[f1*3 + 2] = i1;		
			
			
			setFace((*pIOnum_Face),newPT_ID,v2,i);
			setFace(((*pIOnum_Face)+1),newPT_ID,i1,v2);
			pIOFace[(*pIOnum_Face)*3 + 0] = newPT_ID;
			pIOFace[(*pIOnum_Face)*3 + 1] = v2;
			pIOFace[(*pIOnum_Face)*3 + 2] = i;			
			pIOFace[((*pIOnum_Face)+1)*3 + 0] = newPT_ID;
			pIOFace[((*pIOnum_Face)+1)*3 + 1] = i1;
			pIOFace[((*pIOnum_Face)+1)*3 + 2] = v2;
				

			//find node in cut-path that share face f0 as half edge 
			
			IDCutHedge *nowcut = pCutHedgeH;
			while (next(nowcut) != pCutHedgeT)
			{
				nowcut = next(nowcut);
				if (nowcut->FaceID == f0 && nowcut->NewIDSpilt == v2)
				{		
					nowcut->FaceID = (*pIOnum_Face);
					break;
				}

				if (nowcut->FaceID == f0 && nowcut->NewIDSpilt == i && next(nowcut)->NewIDSpilt == v2)
				{
					next(nowcut)->FaceID = (*pIOnum_Face);
					break;
				}
			}
			
			(*pIOnum_Face)+=2;
			numberF+=2;	
			boundary[newPT_ID] = 0;
		
			//if (count == 3)
			//	return true;
		}
	}
	return change;
}

void	MyParameterization::SortV()
{
#pragma omp parallel for
	for(int i=0;i<numberV;i++)
	{
		if (boundary[i] == 0 )
		{
			if (VHead[i]->next != VTail[i])
			{
				IDList *search = nextN(nextN(VHead[i]));
				IDList *now = search;				
				while (nextN(now) != VTail[i])
				{
					now = nextN(now);
					if (now->ID == search->ID)
					{
						now->back->next = now->next->next;
						now->next->next->back = now->back;

						search->next->back = now->next;
						now->next->next = search->next;
						search->next = now;
						now->back = search;


						search = search->next->next;
						now = search->back;
					}
					else if (now->next->ID == search->ID)
					{
						now->back->next = now->next->next;
						now->next->next->back = now->back;

						//swap now & now->next
						now->back = now->next;
						now->next->next = now;

						search->next->back = now;
						now->next = search->next;
						search->next = now->back;
						search->next->back = search;

						search = search->next->next;
						now = search->back;
					}
					now = next(now);
				}

				//make start node of V same as I
				while (IHead[i]->next->ID != VHead[i]->next->ID)
				{
					IDList *n1 = VHead[i]->next;
					IDList *n2 = VHead[i]->next->next;
					
					VHead[i]->next = n2->next;
					VHead[i]->next->back = VHead[i];
					
					
					n1->back = VTail[i]->back;
					n1->back->next = n1;
					n2->next = VTail[i];
					VTail[i]->back = n2;
				}
				
			}
		}
		else if (boundary[i] == 1 )
		{
			//boundary node
			//find neighbor boundary node
			IDList *now = IHead[i];
			IDList *search = NULL;
			int startID = -1;
			while (nextN(now) != ITail[i])
			{
				now = nextN(now);
				if (boundary[now->ID] == 1)
				{
					startID = now->ID;
					break;
				}
			}
			if (startID < 0)
				throw; //should not happen

			//move boundary VNode to head node
			now = VHead[i];			
			while (next(now) != VTail[i])
			{
				now = nextN(now);
				if (now->ID == startID)
				{
					now->back->next = now->next->next;
					now->next->next->back = now->back;


					now->back = VHead[i];
					now->next->next = VHead[i]->next;
					now->next->next->back = now->next;
					VHead[i]->next = now;
					break;
				}
				now = nextN(now);
				if (now->ID == startID)
				{
					now->back->back->next = now->next;
					now->next->back = now->back->back;

					now->next = now->back;
					now->next->back = now;

					now->back = VHead[i];
					now->next->next = VHead[i]->next;
					now->next->next->back = now->next;
					VHead[i]->next = now;
					break;
				}
			}

			//sort
			search = nextN(nextN(VHead[i]));
			now = search;
			while (nextN(now) != VTail[i])
			{
				now = nextN(now);
				if (now->ID == search->ID)
				{
					now->back->next = now->next->next;
					now->next->next->back = now->back;

					search->next->back = now->next;
					now->next->next = search->next;
					search->next = now;
					now->back = search;


					search = search->next->next;
					now = search->back;
				}
				else if (now->next->ID == search->ID)
				{
					now->back->next = now->next->next;
					now->next->next->back = now->back;

					//swap now & now->next
					now->back = now->next;
					now->next->next = now;

					search->next->back = now;
					now->next = search->next;
					search->next = now->back;
					search->next->back = search;

					search = search->next->next;
					now = search->back;
				}
				now = nextN(now);
			}

			now =now;
		}
		
	}
}

void	MyParameterization::AddCutBorder(						
											IDCutHedge *cutHedgeH,
											IDCutHedge *cutHedgeT,
											IDList *pAddHead,
											IDList *pAddTail,
											int degree,
											int *numDuplicate
										)
{
	IDCutHedge* now = cutHedgeH;
	IDCutHedge* insertPT = NULL;
	int StartAddID = pAddHead->next->ID;
	
	//reset cuthedge structure .... newSpiltID,duplicate 
	//and search for insert point
	while(next(now)!=cutHedgeT)
	{
		now = next(now);
		if (now->NewIDSpilt == StartAddID)
		{
			insertPT = now;
		}
		now->DuplicateDegree = 0;
		//now->twin = NULL;
		now->ID = now->NewIDSpilt;
		now->NewIDSpilt = SPILT_NOT_DETERMINE;
	}
	cutHedgeT->ID = cutHedgeH->next->ID;
	cutHedgeH->ID = cutHedgeT->back->ID;


	if (insertPT != NULL)
	{
		now = insertPT;
		IDCutHedge *newnodeback = (IDCutHedge *)now->back;
		IDList* newCutNode = pAddHead;
		while (newCutNode->next != pAddTail)
		{
			newCutNode = newCutNode->next;
			IDCutHedge::AppendVF(newCutNode->ID,-1,degree,now);
		}

		newCutNode = pAddTail->back;
		while (newCutNode->back != pAddHead->next)
		{
			newCutNode = newCutNode->back;
			IDCutHedge::AppendVF(newCutNode->ID,-1,degree,now);
		}
		//find face of add half edge.  (reverve way: now->back to nowback->next)			
		IDCutHedge *nowNext = NULL;
		IDCutHedge *nowNextNext = NULL;
		while (now->back != newnodeback)
		{
			now = (IDCutHedge *)now->back;
			nowNext = (IDCutHedge *)now->next;
			nowNextNext = (IDCutHedge *)nowNext->next;				
			//find starthalfege;
			int searchHalfEdge = -1;
			for (int i = 0 ; i < 3 ; i++)
			{
				if (Face[nowNext->FaceID][i] != nowNext->ID && Face[nowNext->FaceID][i] != nowNextNext->ID)
				{
					searchHalfEdge = Face[nowNext->FaceID][i];
					break;
				}
			}

			//determine stophalfEdge
			int stopHalfEdge = now->ID;	

			if (searchHalfEdge == stopHalfEdge)  //if start one is same as stop one
			{
				//same face
				now->FaceID = nowNext->FaceID;
			}
			else
			{
				int excludeFaceID = nowNext->FaceID;
				IDList *nowf = FHead[nowNext->ID];
				while(next(nowf)!=FTail[nowNext->ID])
				{
					nowf = next(nowf);
					if (nowf->ID != excludeFaceID)
					{
						for (int i = 0 ; i < 3; i++)
						{
							if (Face[nowf->ID][i] == searchHalfEdge)
							{									
								//find next searchHalfEdge
								if (Face[nowf->ID][(i+1)%3] != nowNext->ID)
								{
									searchHalfEdge = Face[nowf->ID][(i+1)%3];									

								}
								else if (Face[nowf->ID][(i+2)%3] != nowNext->ID)
								{
									searchHalfEdge = Face[nowf->ID][(i+2)%3];									
								}
								else
								{
									printf("ERROR! CANNOT BE THIS CONDITION!\n");
								}

								if (stopHalfEdge == searchHalfEdge)
								{
									//found!
									now->FaceID = nowf->ID;
								}
								excludeFaceID = nowf->ID;
								nowf = FHead[nowNext->ID];
								break;
							}
						}

						if (stopHalfEdge == searchHalfEdge)
						{
							//found!							
							break;
						}
					}
				}////while(next(nowf)!=FTail[nowNext->ID])
			}
		}				
	}
	else
	{
		printf("ERROR!\n");
	}

	//duplicate ID	calculate
	int numDup = 0;
	now = (IDCutHedge *)cutHedgeH;		
	while (next(now) != cutHedgeT)
	{
		now = next(now);
		IDCutHedge *search = now;		
		while(back(search) != cutHedgeH)
		{
			search = back(search);

			if (search->ID == now->ID)
			{
				now->DuplicateDegree = ((IDCutHedge *)search)->DuplicateDegree + 1;
				numDup++;
				//num_boundarySurfacePolarVertex[current_cut_loop]++;
				break;
			}
		}
	}
	*numDuplicate = numDup;

}

void	MyParameterization::CalculateEdgeLength()
{
#pragma omp parallel for
	for (int i=0;i<numberV;i++)
	{
		IDList *now = IHead[i];
		while(now->next!=ITail[i])
		{
			now = now->next;
			now->length =PT->Distance(point[i],point[now->ID]);
		}
		
	} 
}

void	MyParameterization::ChangeVertexIndexInSurfaceFaces(int oldID, int newID, int faceExcludeID, int startHalfEdge,int stopHalfEdge, int *pIOFace)
{
	IDList *now = FHead[oldID];
	int searchHalfEdge  = startHalfEdge;
	if (startHalfEdge == stopHalfEdge)
		return;
	while(next(now)!=FTail[oldID])
	{
		now = next(now);
		if (now->ID != faceExcludeID)
		{
			for (int i = 0 ; i < 3; i++)
			{
				if (Face[now->ID][i] == searchHalfEdge)
				{
					int d[3] = {Face[now->ID][0],Face[now->ID][1],Face[now->ID][2]};					
					//find next searchHalfEdge
					if (Face[now->ID][(i+1)%3] != oldID)
					{
						searchHalfEdge = Face[now->ID][(i+1)%3];
						d[(i+2)%3] = newID;

					}
					else if (Face[now->ID][(i+2)%3] != oldID)
					{
						searchHalfEdge = Face[now->ID][(i+2)%3];
						d[(i+1)%3] = newID;
					}
					else
					{
						printf("ERROR! CANNOT BE THIS CONDITION!\n");
					}
					int recreateFaceID = now->ID;
					DeleteFace(now->ID);

					setFace(recreateFaceID,d[0],d[1],d[2]);
					if (pIOFace)
					{
						pIOFace[(recreateFaceID*3) + 0] = d[0];
						pIOFace[(recreateFaceID*3) + 1] = d[1];
						pIOFace[(recreateFaceID*3) + 2] = d[2];
					}
					if (stopHalfEdge == searchHalfEdge)
						return;
					else
						now = FHead[oldID];
					break;
				}
			}			
		}
	}




}


double	MyParameterization::Parameterize(
											PolarVertex *pIPV,
											int num_PV,
											FILE *logFile)
{
	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;	
    boundarytype=0;	


	//param(logFile);

	MyBoundaryMap();
	//BoundaryMap();
	setPolarMap();
	ParametrizationOptimal(iteNum,PCBCGerror,logFile);
	//ParametrizationSmoothOptimal(iteNum,PCBCGerror,logFile);

	for (int i=0;i<numberV;i++)
	{
		pIPV[i].u = pU[i];
		pIPV[i].v = pV[i];		
	}
	return resultStretch;
}


double	MyParameterization::CircularParameterize(PolarVertex *pIPV,
												 int num_PV,
												 FILE* logFile)
{
	weighttype = 0;
	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;	
    boundarytype=1;	


	//param(logFile);

	MyBoundaryMap();	
	setPolarMap();

	ParametrizationOptimal(iteNum,PCBCGerror,logFile);
	

	for (int i=0;i<numberV;i++)
	{
		pIPV[i].u = pU[i];
		pIPV[i].v = pV[i];		
	}

	getCurrentE();
	return resultStretch;
}


double	MyParameterization::CircularParameterizeOptimalEx(PolarVertex *pIPV,
										int num_PV,
										FILE* logFile)
{
	weighttype = 0;
	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;	
    boundarytype=1;	


	//param(logFile);

	MyBoundaryMap();
	setPolarMap();
	
	
	
	gammaP = 0.75;
	ParametrizationOptimal(iteNum,PCBCGerror,NULL);
	gammaP = 1.0;
	
	

	ParametrizationSmoothOptimal_EX(0.5,0.1,iteNum,PCBCGerror,logFile);

			
	

	for (int i=0;i<numberV;i++)
	{
		pIPV[i].u = pU[i];
		pIPV[i].v = pV[i];		
	}
	return resultStretch;
}
#ifdef INTEL_MKL_VERSION
void	MyParameterization::mkl_Solve()
{
//#define INDEX0
	int nonzero = numberV;
	double *rhs = new double[2*(numberV)];	
	double *solution =  new double[2*(numberV)];
	for(int i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{			
			rhs[i] = 0.0;
			rhs[i+numberV] = 0.0;	//for v
			solution[i] = 0.5;
			solution[i+numberV] = 0.5;
			nonzero += neighborI[i];
		}
		else
		{
			rhs[i] = pU[i];			
			rhs[i+numberV] = pV[i];	//for v
			solution[i] = pU[i];
			solution[i+numberV] = pV[i];
		}
	}	

	int idx = 0;
#ifdef INDEX0
	cpuCompressedMatrixTypeIndex0 spaseMatrixA( numberV*2, numberV,nonzero);
#else
	cpuCompressedMatrixTypeIndex1 spaseMatrixA( numberV*2, numberV,nonzero);
	idx = 1;
#endif
	

	IDList *now = NULL;	
	PolarList *nowp = NULL;
	for(int i=0;i<numberV;i++)
	{
		spaseMatrixA(i+idx,i+idx) = 1.0;
		spaseMatrixA(i+idx +numberV,i+idx) = 1.0;
		
		if(boundary[i]!=1)
		{
			nowp = PHead[i];
      
			while(nextN(nowp)!=PTail[i])
			{
				nowp = nextN(nowp);		
				spaseMatrixA(i+idx,nowp->ID+idx) = -nowp->lambda;
				spaseMatrixA(i+idx + numberV,nowp->ID+idx) = -nowp->lambda;				
			}
		}
		
	}

	MKL_INT *ia = &(spaseMatrixA.index1_data()[0]);
	MKL_INT *ja = &(spaseMatrixA.index2_data()[0]);
	double* eleA  = &(spaseMatrixA.value_data()[0]); 
	
	MKL_INT n = 2*(numberV), rci_request, itercount, expected_itercount = 8;
	MKL_INT ipar[128];
	double *tmp = new double[n*4];
	double dpar[128];
	dcg_init (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
	if (rci_request != 0)
	{
		//fail
		goto exit_function;
	}
		

	ipar[8] = 1;
	ipar[9] = 0;
	dpar[0] = PCBCGerror;

	dcg_check (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
	if (rci_request != 0)
	{
		//failure;
		goto exit_function;
	}


exit_function:	
	delete [] tmp;
	delete [] rhs;
	delete [] solution;

}
#endif
void	MyParameterization::linbcg_Solve()
{
	int i = 0;
	IDList *now = NULL;
	IDList *now2 = NULL;
	PolarList *nowp = NULL;

	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;
	double *UaXY = new double[2*(numberV)+1];
	double *vecb = new double[2*(numberV)+1];
	int nonzero=(numberV);
	for(i=0;i<numberV;i++)
	{
		
		if(boundary[i]!=1)
		{
			nonzero += neighborI[i];			
		}
	}

	int iter=0;
	double linerr=0.0;
	double weight=0.0;
   
	PCBCGSolver *mybcg = new PCBCGSolver(2*nonzero);   
	double *sigsum = new double[numberV];	

	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			mybcg->sa[i+1] = 1.0;
			mybcg->sa[i+1+numberV] = 1.0;  //for v
			vecb[i+1] = 0.0;
			vecb[i+1+numberV] = 0.0;	//for v
		}
		else
		{
			mybcg->sa[i+1] = 1.0;
			vecb[i+1] = pU[i];
			mybcg->sa[i+1+numberV] = 1.0;	//for v
			vecb[i+1+numberV] = pV[i];	//for v
		}
	}
	mybcg->ija[1] = 2*(numberV)+2;
	int dlk=2*(numberV)+1;
   
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			nowp = PHead[i];
      
			while(nextN(nowp)!=PTail[i])
			{
				nowp = nextN(nowp);
				++dlk;
				mybcg->sa[dlk] = -nowp->lambda;
				mybcg->ija[dlk]=nowp->ID+1;
			}
		}
		mybcg->ija[i+1+1]=dlk+1;
	}

	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			nowp = PHead[i];
			while(nextN(nowp)!=PTail[i])
			{
				nowp = nextN(nowp);
				++dlk;
				mybcg->sa[dlk] = -nowp->lambda;
				mybcg->ija[dlk]=nowp->ID+numberV+1;
			}
		}
		mybcg->ija[i+numberV+1+1]=dlk+1;
	}
  
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]==1)
		{
			UaXY[i+1] = pU[i];
			UaXY[i+numberV+1] = pV[i];
		}
		else
		{
			UaXY[i+1] = 0.5;
			UaXY[i+numberV+1] = 0.5;
		}
	}

	mybcg->linbcg(((unsigned long)(2*(numberV))),vecb,UaXY,1,PCBCGerror,(1000+iteNum),&iter,&linerr);
  
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			pU[i] = UaXY[i+1];
			pV[i] = UaXY[i+numberV+1];
		}
	}
	delete [] UaXY;
	delete [] vecb;
   	delete  mybcg;   
	delete [] sigsum;

}

void    MyParameterization::FindMinMaxStretchBoundaryFace(int *minID,int *maxID)
{
	double pV1,pV2,pV3,pU1,pU2,pU3;
	double dE,dG;
	double dsize1;
	dE = 0.0;
	dG = 0.0;
	double dsum = 0;
	double min_stretch = DBL_MAX;
	double max_stretch = 0.0;
	int max_stretch_face = -1;
	int min_stretch_face = -1;

	for(int i=0;i<numberF;i++)
	{
		if(	boundary[Face[i][0]]==1	||
			boundary[Face[i][1]]==1 ||
			boundary[Face[i][2]]==1)
		{      
			pV1 = pV[Face[i][0]];
			pV2 = pV[Face[i][1]];
			pV3 = pV[Face[i][2]];
			pU1 = pU[Face[i][0]];
			pU2 = pU[Face[i][1]];
			pU3 = pU[Face[i][2]];
			
			dsize1 = PT->getParametricA(pV1,pV2,pV3,pU1,pU2,pU3);
      
			PT->setParametricDs(bc[0],point[Face[i][0]],
					point[Face[i][1]],point[Face[i][2]],
					pV1,pV2,pV3,dsize1);
			PT->setParametricDt(bc[1],point[Face[i][0]],
					point[Face[i][1]],point[Face[i][2]],
					pU1,pU2,pU3,dsize1);
			dE = PT->InnerProduct(bc[0],bc[0]);
    
			dG = PT->InnerProduct(bc[1],bc[1]);

			double stretch = (dE+dG);
			if (stretch > max_stretch)
			{
				max_stretch = stretch;
				max_stretch_face = i;
			}

			if (stretch < min_stretch)
			{
				min_stretch = stretch;
				min_stretch_face = i;
			}
		}
	}

	*minID = min_stretch_face;
	*maxID = max_stretch_face;
}
double   MyParameterization::GetExtremaTriangle(int *index , bool *borderFace,int SolverMode)
{
	
	
	int i = 0;


	boundarytype=1;
	MyBoundaryMap();
	setPolarMap(); 	
	
	//setTutteC();
	setFloaterC();
	//setMVCC();
	//setMVCC();
	SortIndexP();


	linbcg_Solve();
	
	/*

	FILE *fpWeight = fopen("stretch_weight_setting.txt","rt");
	double c = 0;			
	double m = 0;

	if (fpWeight)
	{
		printf("FOUND WEIGHT SETTING FILE!\n");
		//read setting from file
		float input;
		fscanf(fpWeight,"c=%f\r\n",&input);
		c = input;
		//fseek(fpWeight, 0L, SEEK_SET );
		fscanf(fpWeight,"m=%f\r\n",&input);
		m = input;
		fclose(fpWeight);

		printf("c = %f\n"
			   "m = %f\n",c,m);
		 
	}
	*/

	
	//double dE,dG;
	//double dsize1;
	//dE = 0.0;
	//dG = 0.0;
	double dsum = 0;
	double max_stretch = 0.0;
	int max_stretch_face = -1;

	

#pragma omp parallel for
	for(i=0;i<numberF;i++)
	{
		//if(boundary[Face[i][0]]!=1&&boundary[Face[i][1]]!=1&&boundary[Face[i][2]]!=1)
		{
			double pV1,pV2,pV3,pU1,pU2,pU3;
			pV1 = pV[Face[i][0]];
			pV2 = pV[Face[i][1]];
			pV3 = pV[Face[i][2]];
			pU1 = pU[Face[i][0]];
			pU2 = pU[Face[i][1]];
			pU3 = pU[Face[i][2]];
			/*
			double centerU = (pU1+pU2+pU3)/3.0;
			double centerV = (pV1+pV2+pV3)/3.0;
			double distance = sqrt(centerU*centerU + centerV*centerV);
			*/
			double dsize1 = PT->getParametricA(pV1,pV2,pV3,pU1,pU2,pU3);
    
   
			Point3d _bc[2];
			PT->setParametricDs(&_bc[0],point[Face[i][0]],
					point[Face[i][1]],point[Face[i][2]],
					pV1,pV2,pV3,dsize1);
			PT->setParametricDt(&_bc[1],point[Face[i][0]],
					point[Face[i][1]],point[Face[i][2]],
					pU1,pU2,pU3,dsize1);
			double dE = PT->InnerProduct(&_bc[0],&_bc[0]);    
			double dG = PT->InnerProduct(&_bc[1],&_bc[1]);
			double stretch = (dE+dG);

			#pragma omp critical
			{
				dsum +=  areaMap3D[i]*0.5*(stretch);
				if (stretch > max_stretch)
				{
					max_stretch = stretch;
					max_stretch_face = i;
				}
			}

			
		}
	}

	if (max_stretch_face < 0)
	{		
		return -1;
	}
	*index = max_stretch_face;	
	*borderFace = boundary[Face[max_stretch_face][0]] || boundary[Face[max_stretch_face][1]] || boundary[Face[max_stretch_face][2]];
	//boundarysigma=0;
	double E =constsumarea3D*sqrt(dsum/sumarea3D);
	return E;
}


double	MyParameterization::GetExtremaVertex(int *index , bool *borderVertex,int SolverMode)
{
	int i = 0;
	boundarytype=1;
	MyBoundaryMap();
	setPolarMap(); 
	
	
	
	
	//setTutteC();
	setFloaterC();
	//setMVCC();
	
	SortIndexP();
	linbcg_Solve();
	
	
	


	FILE *fpWeight = fopen("stretch_weight_setting.txt","rt");
	double c = 0;			
	double m = 0;

	if (fpWeight)
	{
		printf("FOUND WEIGHT SETTING FILE!\n");
		//read setting from file
		float input;
		fscanf(fpWeight,"c=%f\r\n",&input);
		c = input;
		//fseek(fpWeight, 0L, SEEK_SET );
		fscanf(fpWeight,"m=%f\r\n",&input);
		m = input;
		fclose(fpWeight);

		printf("c = %f\n"
			   "m = %f\n",c,m);
		 
	}
	double max_stretch = 0.0;
	int max_stretch_face = -1;
	
	IDList *now=NULL;
	double varphi,ddv,dsize1,sumarea;
	double dddhval=0.0;
	double localsum=0.0;
	for(i=0;i<numberF;i++)
	{
    
		dsize1 = PT->getParametricA(pV[Face[i][0]],
					pV[Face[i][1]],
					pV[Face[i][2]],
					pU[Face[i][0]],
					pU[Face[i][1]],
					pU[Face[i][2]]);
		PT->setParametricDs(bc[0],point[Face[i][0]],
				point[Face[i][1]],point[Face[i][2]],
				pV[Face[i][0]],pV[Face[i][1]],
				pV[Face[i][2]],dsize1);
		PT->setParametricDt(bc[1],point[Face[i][0]],
				point[Face[i][1]],point[Face[i][2]],
				pU[Face[i][0]],pU[Face[i][1]],
				pU[Face[i][2]],dsize1);
    
		E[i] = PT->InnerProduct(bc[0],bc[0]);
		G[i] = PT->InnerProduct(bc[1],bc[1]);

		double stretch = E[i]+G[i];		
		if (stretch > max_stretch)
		{
			max_stretch = stretch;
			max_stretch_face = i;
		}
	}

	int maxStretchIndex = 0;
	double local_sigma[3] = {0};
	for(i=0;i<3;i++)
	{         
		int index = Face[max_stretch_face][i];
		local_sigma[i]=0.0;
		now = FHead[index];
		varphi=0.0;
		localsum=0.0;    
		while(next(now)!=FTail[index])
		{
			now = next(now);
			varphi += (areaMap3D[now->ID]*(0.5*(E[now->ID]+G[now->ID])));
			localsum += (areaMap3D[now->ID]);      
		}   
		local_sigma[i] = sqrt((varphi/localsum));    

		if (local_sigma[i] > local_sigma[maxStretchIndex])
			maxStretchIndex = i;
      
	}

	*index = Face[max_stretch_face][maxStretchIndex];
	if ( boundary[*index]  == 1)
		*borderVertex = true;
	else
		*borderVertex = false;

	double E = getCurrentE();
	//boundarysigma=1;
	return E;
}


/*
void	MyParameterization::Parameterize(
											PolarVertex *pIPV,
											int num_PV,
											int	*pInputFace,
											int num_Face
										)
{
	memoryallocate( num_PV, num_Face);
	for(int i=0;i<num_PV;i++)
	{		
		setPoint(i,pIPV[i].p_vertex->x,pIPV[i].p_vertex->y,pIPV[i].p_vertex->z);
	}
	
	iteNum = ((num_PV/20000) + 1)*2000;
	for(int i=0;i<num_Face;i++)
	{		
		setFace(i,pInputFace[(i*3)+0],pInputFace[(i*3)+1],pInputFace[(i*3)+2]);
		IDtool->AppendVFSort(i,FHead[Face[i][0]],FTail[Face[i][0]]);
		IDtool->AppendVFSort(i,FHead[Face[i][1]],FTail[Face[i][1]]);
		IDtool->AppendVFSort(i,FHead[Face[i][2]],FTail[Face[i][2]]);
	}   
	// feature analysis 
    
	SetBoundaryLines();
	setAreaMap3D();

	param();

	//output
	for (int i=0;i<numberV;i++)
	{
		pIPV[i].u = pU[i];
		pIPV[i].v = pV[i];		
	}
	
}
*/

//if changed or new data.... return true;
bool    MyParameterization::InsertDijkstra(int ID, double length,IDList *pHead,IDList *pTail)
{
	IDList *search = IDtool->GetI(ID,pHead,pTail);
	if (search)
	{
		if (search->length > length)
		{
			search->length = length;			
			IDList *dummy = search->back;
			if(dummy->back!=pHead)
			{
				while(search->length < dummy->length)
				{
					dummy = dummy->back;
					if(dummy==pHead)
						break;
				}
			}

			if (dummy != search->back)
			{
				search->back->next = search->next;
				search->next->back = search->back;
				IDList *dummynext=dummy->next;      
				search->back = dummy;
				dummy->next = search;
				search->next = dummynext;
				dummynext->back =search;
			}
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		//not found
		IDList *now = new IDList(ID);
		now->length = length;
		IDList *dummy = pHead->next;
		if(pHead->next!=pTail)
		{
			while(now->length > dummy->length)
			{
				dummy = dummy->next;
				if(dummy==pTail)
					break;
			}
		}
		IDList *dummyback=dummy->back;
      
		now->next = dummy;
		dummy->back = now;
		now->back = dummyback;
		dummyback->next =now;
		return true;	
	}
	
}


int	MyParameterization::FindShortestPathToBorderFromPoint( int startVertexID,IDList *pHead,IDList *pTail,double *length)
{
	int *previous = new int [numberV];
	memset(previous,0xFF,sizeof(int)*numberV);

	bool *check = new bool [numberV];
	memset(check,0,sizeof(bool)*numberV);

	IDList *pDijHead = new IDList();
	IDList *pDijTail = new IDList();

	pDijHead->next = pDijTail;
	pDijTail->back = pDijHead;

	InsertDijkstra(startVertexID,0.0,pDijHead,pDijTail);
	IDList *smallest = NULL;
	while(1)
	{
		smallest = pDijHead->next;
		if (smallest == pDijTail)
			break; //not found shortest path??

		if (boundary[smallest->ID] == 1)
			break;

		pDijHead->next = smallest->next;
		smallest->next->back = pDijHead;
		check[smallest->ID] = true;
		IDList *now = IHead[smallest->ID];
		while (now->next != ITail[smallest->ID])
		{
			now = now->next;			
			
			if (!check[now->ID])
			{
				if (InsertDijkstra(now->ID,now->length + smallest->length,pDijHead,pDijTail))
				{
					previous[now->ID] = smallest->ID;
				}
			}
		}
		delete smallest;
	}

	int numNode = 0;
	//create output
	int currentID = smallest->ID;

	if (currentID == startVertexID)
	{
		if (length)
			*length = 0;
		IDtool->CleanNeighbor(pDijHead,pDijTail);
		delete [] previous;
		delete [] check;
		return numNode;
	}
	if (length)
		*length = smallest->length;
	while (currentID != startVertexID)
	{
		IDtool->AppendVF(currentID,pTail);		 
		currentID = previous[currentID];
		numNode++;
	}
	IDtool->AppendVF(startVertexID,pTail);	
	numNode++;

	IDtool->CleanNeighbor(pDijHead,pDijTail);
	delete [] previous;
	delete [] check;
	return numNode;

}



int	MyParameterization::FindShortestPathToBorderFromFace( int startFaceID,IDList *pHead,IDList *pTail)
{
	int *previous = new int [numberV];
	memset(previous,0xFF,sizeof(int)*numberV);

	bool *check = new bool [numberV];
	memset(check,0,sizeof(bool)*numberV);

	IDList *pDijHead = new IDList();
	IDList *pDijTail = new IDList();

	pDijHead->next = pDijTail;
	pDijTail->back = pDijHead;


	InsertDijkstra(Face[startFaceID][0],0.0,pDijHead,pDijTail);
	InsertDijkstra(Face[startFaceID][1],0.0,pDijHead,pDijTail);
	InsertDijkstra(Face[startFaceID][2],0.0,pDijHead,pDijTail);
	IDList *smallest = NULL;
	while(1)
	{
		smallest = pDijHead->next;
		if (smallest == pDijTail)
			break; //not found shortest path??

		if (boundary[smallest->ID] == 1)
			break;

		pDijHead->next = smallest->next;
		smallest->next->back = pDijHead;
		check[smallest->ID] = true;
		IDList *now = IHead[smallest->ID];
		while (now->next != ITail[smallest->ID])
		{
			now = now->next;			
			
			if (!check[now->ID])
			{
				if (InsertDijkstra(now->ID,now->length + smallest->length,pDijHead,pDijTail))
				{
					previous[now->ID] = smallest->ID;
				}
			}
		}
		delete smallest;
	}

	int numNode = 0;
	//create output
	int currentID = smallest->ID;

	if (currentID == Face[startFaceID][0] ||currentID == Face[startFaceID][1] ||currentID == Face[startFaceID][2] )
	{
		IDtool->CleanNeighbor(pDijHead,pDijTail);
		delete [] previous;
		delete [] check;
		return numNode;
	}

	while (currentID != Face[startFaceID][0] && currentID != Face[startFaceID][1] && currentID != Face[startFaceID][2])
	{
		IDtool->AppendVF(currentID,pTail);		 
		currentID = previous[currentID];
		numNode++;
	}
	IDtool->AppendVF(currentID,pTail);	
	numNode++;

	IDtool->CleanNeighbor(pDijHead,pDijTail);
	delete [] previous;
	delete [] check;

	return numNode;

}

void MyParameterization::setTutteC()
{
	int i;	
	IDList *now=NULL;
	PolarList *nowp=NULL;
	PolarList *nowp2=NULL;	
	double dweight=0.0;

	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			dweight=0.0;
			nowp = PHead[i];
			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);
				dweight += 1.0;
				nowp->lambda += 1.0;				
			}

			if(dweight!=0.0)
				nowp = PHead[i];

			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);
				nowp->lambda /= dweight;
			}
		}
	}
}
void    MyParameterization::CalBorderPath(IDList *BpointH,IDList *BpointT,double *length,int *numPoint)
{
	int i=0,j=0,k=0;
	IDList *now=NULL;
	IDList *BedgeH = new IDList();
	IDList *BedgeT = new IDList();
	BedgeH->next = BedgeT;
	BedgeT->back = BedgeH;
	int cnt=0;
    
	for(i=0;i<numberF;i++)
	{
		if(	boundary[Face[i][0]]==1&&
			boundary[Face[i][1]]==1&&
			boundary[Face[i][2]]==1)
		{
	
			IDtool->AppendVF(Face[i][0],BedgeT);
			IDtool->AppendVF(Face[i][1],BedgeT);
			IDtool->AppendVF(Face[i][1],BedgeT);
			IDtool->AppendVF(Face[i][2],BedgeT);
			IDtool->AppendVF(Face[i][2],BedgeT);
			IDtool->AppendVF(Face[i][0],BedgeT);
		}
		else
		{
			if(	(boundary[Face[i][0]]==1&&boundary[Face[i][1]]==1)||
				(boundary[Face[i][1]]==1&&boundary[Face[i][2]]==1)||
				(boundary[Face[i][2]]==1&&boundary[Face[i][0]]==1))
			{
				if(boundary[Face[i][0]]==1&&boundary[Face[i][1]]==1)
				{
					IDtool->AppendVF(Face[i][0],BedgeT);
					IDtool->AppendVF(Face[i][1],BedgeT);	    
				}
				else if(boundary[Face[i][1]]==1&&boundary[Face[i][2]]==1)
				{
					IDtool->AppendVF(Face[i][1],BedgeT);
					IDtool->AppendVF(Face[i][2],BedgeT);	    
				}
				else if(boundary[Face[i][2]]==1&&boundary[Face[i][0]]==1)
				{
					IDtool->AppendVF(Face[i][2],BedgeT);
					IDtool->AppendVF(Face[i][0],BedgeT);	    
				}
			}
		}
	}
    

	//filter of real border edge from a pair of vertex that indicated as border point
	int jumpcheck=0;
	now = BedgeH;
	while(next(now)!=BedgeT)
	{
		now = next(now);
		jumpcheck=0;
		if(boundary[now->ID]== 1&&boundary[next(now)->ID]==1)
		{
			int checkbedge = 0;
			IDList *now2 = VHead[now->ID];
			while(next(now2)!=VTail[now->ID])
			{
				now2 = next(now2);
				if(	now2->ID==next(now)->ID||
					next(now2)->ID==next(now)->ID)
				{
					checkbedge++;
				}
				now2 = next(now2);
			}

			if(checkbedge>=2)
			{	
				//remove this edge (pair of vertice) since this edge is not really border.

				IDList *dummyn = next(now)->next;
				IDList *dummyb = now->back;
				dummyb->next = dummyn;
				dummyn->back = dummyb;
	  
				IDList *ddn = next(now);
				IDList *ddnn = now;
				delete ddn;
				delete ddnn;
	  
				now = dummyb;
				jumpcheck=1;
			}
		}

		if(jumpcheck!=1)
			now = next(now);
	}

	double tlength=0.0;
	now = BedgeH;
	while(next(now)!=BedgeT)
	{
		now = next(now);
		tlength += PT->Distance(point[now->ID],point[next(now)->ID]);
		now = next(now);
      
	}
	
	//tlength is total length of border now

	int startID=0;
	int *checkbin = new int[numberV];
	for(i=0;i<numberV;i++)
		checkbin[i]=-1;

	//select start point of border
	if(pickID==-1)
	{
		for(i=0;i<numberV;i++)
		{
			if(boundary[i]==1)
			{
				pickID=i;
				break;
			}
		}
	}
    
	//printf("start id = %d\n",pickID);
	startID=pickID;
	checkbin[pickID]=0;
    
	IDtool->AppendVF(pickID,BpointT);
	k=1;
  
	if (numPoint)
		*numPoint = 0;
	while(1)
	{
		now = BedgeH;
    
		while(next(now)!=BedgeT)
		{
			now = next(now);
			if(startID==now->ID||startID==next(now)->ID)
			{
				if(startID==now->ID&&checkbin[next(now)->ID]==-1)
				{
					k++;
					checkbin[next(now)->ID]=0;
					IDtool->AppendVF(next(now)->ID,BpointT);
					*numPoint += 1;
					startID=next(now)->ID;
					break;
				}
				else if(startID==next(now)->ID&&checkbin[now->ID]==-1)
				{
					k++;
					checkbin[now->ID]=0;
					IDtool->AppendVF(now->ID,BpointT);
					*numPoint += 1;
					startID=now->ID;
					break;
				}
	
			}
			now = next(now);
		}
    
		if(k>=numboundary)
			break;
	}
  
	IDtool->CleanNeighbor(BedgeH,BedgeT);  //de-allocate BedgeH&T because border info now is at BpointH&T
	delete [] checkbin;
	*length = tlength;
}

void MyParameterization::MyBoundaryMap()
{
	int i=0,j=0,k=0;
  
	if (pU)
		delete [] pU;
	if (pV)
		delete [] pV;

	pU = new double[numberV];
	pV = new double[numberV];
	memset(pU,0x00,sizeof(double)*numberV);
	memset(pV,0x00,sizeof(double)*numberV);

    
	IDList *now=NULL;
	IDList *BedgeH = new IDList();
	IDList *BedgeT = new IDList();
	BedgeH->next = BedgeT;
	BedgeT->back = BedgeH;
	int cnt=0;
 #pragma omp parallel for			   
	for(i=0;i<numberF;i++)
	{
		if(	boundary[Face[i][0]]==1&&
			boundary[Face[i][1]]==1&&
			boundary[Face[i][2]]==1)
		{
			#pragma omp critical
			{
				IDtool->AppendVF(Face[i][0],BedgeT);
				IDtool->AppendVF(Face[i][1],BedgeT);
				IDtool->AppendVF(Face[i][1],BedgeT);
				IDtool->AppendVF(Face[i][2],BedgeT);
				IDtool->AppendVF(Face[i][2],BedgeT);
				IDtool->AppendVF(Face[i][0],BedgeT);
			}
		}
		else
		{
			if(	(boundary[Face[i][0]]==1&&boundary[Face[i][1]]==1)||
				(boundary[Face[i][1]]==1&&boundary[Face[i][2]]==1)||
				(boundary[Face[i][2]]==1&&boundary[Face[i][0]]==1))
			{
				
				
				if(boundary[Face[i][0]]==1&&boundary[Face[i][1]]==1)
				{
					#pragma omp critical
					{

						IDtool->AppendVF(Face[i][0],BedgeT);
						IDtool->AppendVF(Face[i][1],BedgeT);	    
					}
				}
				else if(boundary[Face[i][1]]==1&&boundary[Face[i][2]]==1)
				{
					#pragma omp critical
					{
					IDtool->AppendVF(Face[i][1],BedgeT);
					IDtool->AppendVF(Face[i][2],BedgeT);	    
					}
				}
				else if(boundary[Face[i][2]]==1&&boundary[Face[i][0]]==1)
				{
					#pragma omp critical
					{
					IDtool->AppendVF(Face[i][2],BedgeT);
					IDtool->AppendVF(Face[i][0],BedgeT);
					}
				}
				
			}
		}
	}
    

	//filter of real border edge from a pair of vertex that indicated as border point
	int jumpcheck=0;
	now = BedgeH;
	while(next(now)!=BedgeT)
	{
		now = nextN(now);
		jumpcheck=0;
		if(boundary[now->ID]== 1&&boundary[nextN(now)->ID]==1)
		{
			int checkbedge = 0;
			IDList *now2 = VHead[now->ID];
			while(nextN(now2)!=VTail[now->ID])
			{
				now2 = nextN(now2);
				if(	now2->ID==nextN(now)->ID||
					nextN(now2)->ID==nextN(now)->ID)
				{
					checkbedge++;
				}
				now2 = nextN(now2);
			}

			if(checkbedge>=2)
			{	
				//remove this edge (pair of vertice) since this edge is not really border.

				IDList *dummyn = nextN(now)->next;
				IDList *dummyb = now->back;
				dummyb->next = dummyn;
				dummyn->back = dummyb;
	  
				IDList *ddn = nextN(now);
				IDList *ddnn = now;
				delete ddn;
				delete ddnn;
	  
				now = dummyb;
				jumpcheck=1;
			}
		}

		if(jumpcheck!=1)
			now = nextN(now);
	}
    

    
    
	IDList *BpointH = new IDList();
	IDList *BpointT = new IDList();
    
	BpointH->next = BpointT;
	BpointT->back = BpointH;
	double tlength=0.0;
	now = BedgeH;
	while(nextN(now)!=BedgeT)
	{
		now = nextN(now);
		tlength += PT->Distance(point[now->ID],point[next(now)->ID]);
		now = nextN(now);
      
	}
	
	//tlength is total length of border now

	int startID=0;
	int *checkbin = new int[numberV];
	for(i=0;i<numberV;i++)
		checkbin[i]=-1;

	//select start point of border
	if(pickID==-1)
	{
		for(i=0;i<numberV;i++)
		{
			if(boundary[i]==1)
			{
				pickID=i;
				break;
			}
		}
	}
    
	//printf("start id = %d\n",pickID);
	startID=pickID;
	checkbin[pickID]=0;
    
	IDtool->AppendVF(pickID,BpointT);
	k=1;
  
  
	while(1)
	{
		now = BedgeH;
    
		while(nextN(now)!=BedgeT)
		{
			now = nextN(now);
			if(startID==now->ID||startID==nextN(now)->ID)
			{
				if(startID==now->ID&&checkbin[next(now)->ID]==-1)
				{
					k++;
					checkbin[nextN(now)->ID]=0;
					IDtool->AppendVF(nextN(now)->ID,BpointT);
					startID=nextN(now)->ID;
					break;
				}
				else if(startID==nextN(now)->ID&&checkbin[now->ID]==-1)
				{
					k++;
					checkbin[now->ID]=0;
					IDtool->AppendVF(now->ID,BpointT);
					startID=now->ID;
					break;
				}
	
			}
			now = nextN(now);
		}
    
		if(k>=numboundary)
			break;
	}
  
	IDtool->CleanNeighbor(BedgeH,BedgeT);  //de-allocate BedgeH&T because border info now is at BpointH&T
	
  
	double clen=0.0;
	
	i=0;
  
  
	/* Unit Circle */
	if(boundarytype==1)
	{
		double length = 0;
		double unitArc = (2*PI)/tlength;
		now = next(BpointH);
		pU[now->ID] = 0.5*cos(0.0);
		pV[now->ID] = 0.5*sin(0.0);
		while(next(now)!=BpointT)
		{
			now = next(now);
			length += PT->Distance(point[now->ID],point[back(now)->ID]);
			//pU[now->ID] = 0.5*cos(2.0*PI*(((double)(i))/((double)(numboundary))));
			//pV[now->ID] = 0.5*sin(2.0*PI*(((double)(i))/((double)(numboundary))));
			pU[now->ID] = 0.5*cos(unitArc*length);
			pV[now->ID] = 0.5*sin(unitArc*length);
			i++;
		}
	}
	else
	{    
		/* Unit Square */
		tlength *= 0.25;
		now = BpointH;
		BpointT->ID = BpointH->next->ID; // for cal length;
		int state=0;
		clen=0.0;
		double this_side_length[4] = {0};		
		
		int avg = numboundary/4;
		if (avg == 0) 
			avg = 1;
		/*
		switch(numboundary%4)
		{
			case 0:
				this_side_num_edge[0] = avg;
				this_side_num_edge[1] = avg;
				this_side_num_edge[2] = avg;
				this_side_num_edge[3] = avg;
				break;
			case 1:
				this_side_num_edge[0] = avg+1;
				this_side_num_edge[1] = avg;
				this_side_num_edge[2] = avg;
				this_side_num_edge[3] = avg;
				break;
			case 2:
				this_side_num_edge[0] = avg+1;
				this_side_num_edge[1] = avg;
				this_side_num_edge[2] = avg+1;
				this_side_num_edge[3] = avg;
				break;
			case 3:
				this_side_num_edge[0] = avg+1;
				this_side_num_edge[1] = avg+1;
				this_side_num_edge[2] = avg+1;
				this_side_num_edge[3] = avg;
				break;

		}
		*/
		double sum_length = 0;
		double edge_length = 0;
		double loop = 0;
		int min_error_edgediff = numboundary*numboundary;
		int min_error_id = -1;
		IDList *startPoint = BpointH;
		
		
		while (loop < tlength)
		{
			state = 0;
			sum_length = 0;
			double temp_length[4] ={0};
			int this_side_num_edge[4] = {0};
			
			startPoint = next(startPoint);
			IDList *now = startPoint;
			do
			{				
				edge_length = PT->Distance(point[now->ID],point[next(now)->ID]);			
				if (sum_length + edge_length>= tlength && state < 3)
				{
					if ((sum_length + edge_length)/tlength >= 1.1&& this_side_num_edge[state] > 0)					
					{
						temp_length[state] = sum_length;
						this_side_num_edge[state+1]++;
						sum_length = edge_length;
					}
					else
					{
						temp_length[state] = sum_length+edge_length;
						this_side_num_edge[state]++;
						sum_length = 0;
					}				
					state++;
				}
				else
				{
					this_side_num_edge[state]++;
					sum_length += edge_length;
				}
				now = next(now);
				if (now == BpointT)
					now = BpointH->next;
			}
			while (now != startPoint);			
			temp_length[state] = sum_length;

			//cal error
			int thisError = 0;
			for (int e = 0 ; e < 4 ; e++)
			{
				if (this_side_num_edge[e] > 0)
					thisError += (avg-this_side_num_edge[e])*(avg-this_side_num_edge[e]);
			}

			if (thisError < min_error_edgediff)
			{
				min_error_id = startPoint->ID;
				min_error_edgediff= thisError;
				this_side_length[0] = temp_length[0];
				this_side_length[1] = temp_length[1];
				this_side_length[2] = temp_length[2];
				this_side_length[3] = temp_length[3];
			}
			loop += PT->Distance(point[startPoint->ID],point[next(startPoint)->ID]);
		}
		
		//reconstruct BpointH , BpointT according to min_error_id		
		while (BpointH->next->ID != min_error_id)
		{
			IDList *move_node = BpointH->next;
			BpointH->next = move_node->next;
			move_node->next->back = BpointH; 
			
			move_node->back = BpointT->back;
			BpointT->back->next = move_node;			
			BpointT->back = move_node;
			move_node->next = BpointT;			
		}

		printf("TEST POINT CORNER ID: %d ",BpointH->next->ID);

		state = 0;
		sum_length  = 0;
		now = BpointH;
		BpointT->ID = BpointH->next->ID; // for cal length;
		while(next(now)!=BpointT)
		{
			now = next(now);
			switch(state)
			{
				case 0:
					if (clen < 1.0)
					{
						pU[now->ID] = clen; 
						pV[now->ID] = 0.0;
					}
					else
					{
						pU[now->ID] = 1.0;
						pV[now->ID] = 0.0;
						sum_length = 0.0;
						state=1;
						printf("%d ",now->ID);
						
					}
					break;
				case 1:
					if (clen < 1.0)
					{
						pU[now->ID] = 1.0; 
						pV[now->ID] = clen;
					}
					else
					{
						pU[now->ID] = 1.0;
						pV[now->ID] = 1.0;
						sum_length = 0.0;
						state=2;
						printf("%d ",now->ID);

					}
					break;
				case 2:
					if (clen < 1.0)
					{
						pU[now->ID] = 1.0-clen; 
						pV[now->ID] = 1.0;
					}
					else
					{
						pU[now->ID] = 0.0;
						pV[now->ID] = 1.0;
						sum_length = 0.0;
						state=3;
						printf("%d\n",now->ID);

					}
					break;
				case 3:
					if (clen < 1.0)
					{
						pU[now->ID] = 0.0; 
						pV[now->ID] = 1.0-clen;
					}
					else
					{
						//should never enter this one;
						state=4;
					}
				default:
					break;

			}		
			sum_length += PT->Distance(point[now->ID],point[next(now)->ID]);
			clen = sum_length/this_side_length[state];			
			//sum_length += 1.0;
			//clen = sum_length/this_side_num_edge[state];			
		}		
		
	}
double init_uncontraintsU;
double init_uncontraintsV;

if(boundarytype==1)
{
	init_uncontraintsU = 0.5;
	init_uncontraintsV = 0.5;
}
else
{
	init_uncontraintsU = 0.0;
	init_uncontraintsV = 0.0;
}

//#pragma omp parallel for
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]==1)
		{            
		
		}
		else
		{
			pU[i] = init_uncontraintsU;
			pV[i] = init_uncontraintsV;
		}
	}
  
	IDtool->CleanNeighbor(BpointH,BpointT);
	delete [] checkbin;  

	if(boundarytype==0||boundarytype==2){
    constsumarea3D = sqrt(1.0/sumarea3D);
  }else if(boundarytype==1){
    constsumarea3D = sqrt(((0.5*0.5*PI)/sumarea3D));
  }

}

double	MyParameterization::ParamWithFaceNormalStretch(	PolarVertex *pIPV,
														int num_PV,
														FILE* logFile)
{
	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;	
    boundarytype=0;
	MyBoundaryMap();
	//BoundaryMap();
	setPolarMap();
	ParametrizationOptimal(iteNum,PCBCGerror,logFile);
	/*
	for (int i=0;i<numberV;i++)
	{
		pIPV[i].u = pU[i];
		pIPV[i].v = pV[i];
	}	
	return 0;
	*/

	double *UaXY = new double[2*(numberV)+1];
    double *vecb = new double[2*(numberV)+1];  
	int i;
	IDList *now;
	IDList *now2;
	PolarList *nowp;
	level=0;  
	int nonzero=(numberV);
	for(i=0;i<numberV;i++)
	{
		vecb[i+1]=0.0;
		if(boundary[i]!=1)
		{
			now = IHead[i];
			while(next(now)!=ITail[i])
			{
				now = next(now);
				nonzero++;
			}
		}
	}

	int iter=0;
	double linerr=0.0;
	double weight=0.0;
  
	PCBCGSolver *mybcg = new PCBCGSolver(2*nonzero);
	double *sigsum = new double[numberV];
	memset(sigsum,0x00,sizeof(double)*numberV);

	PointTool pt;
	setSigmaZero();
	for(i=0;i<numberV;i++)
	{      
		if(boundary[i]!=1)
		{
			IDList *nowv = VHead[i];
			double point_normal[3] = {0};
			int numSurroundFace = neighborF[i];
			Point3d *faceNormal = new Point3d[numSurroundFace];
			Point3d pointNormal;
			Point3d v1,v2,N;
			int j = 0;
			while (next(nowv) != VTail[i])
			{
				nowv = next(nowv);
				pt.makeVector(&v1,point[i],point[nowv->ID]);
				pt.makeVector(&v2,point[i],point[nowv->next->ID]);
				pt.CrossVector(&N,&v1,&v2);
				faceNormal[j].x =  N.x;faceNormal[j].y =  N.y;faceNormal[j].z =  N.z;
				j++;
				pointNormal.x += N.x;pointNormal.y += N.y;pointNormal.z += N.z;
				nowv = next(nowv);
			}
			pt.Normalize3D(&pointNormal);
					
		
			nowp = PHead[i];
			sigsum[i]=0.0;
			j = numSurroundFace-1;
			while (next(nowp) != PTail[i])
			{
				nowp = next(nowp);
				// faceNormal[j%numSurroundFace]; //v1
				// faceNormal[(j+1)%numSurroundFace];	 //v2
				N.x = faceNormal[j%numSurroundFace].x + faceNormal[(j+1)%numSurroundFace].x;
				N.y = faceNormal[j%numSurroundFace].y + faceNormal[(j+1)%numSurroundFace].y;
				N.z = faceNormal[j%numSurroundFace].z + faceNormal[(j+1)%numSurroundFace].z;
				pt.Normalize3D(&N);
				double cosAngle = pt.InnerProduct(&pointNormal,&N);
				if (cosAngle > 1.0)
					cosAngle = 1.0;
				if (cosAngle < -1.0)
					cosAngle = -1.0;
				double Angle = acos(cosAngle);
				nowp->lambda = nowp->old_lambda;
				
				
				/*
				if (sig < 1)
				{
					sig /= (2 - cosAngle);
				}
				else
				{
					sig *= (2 - cosAngle);
				}
				*/
				//nowp->lambda /= sig;
				
				nowp->lambda *= pow((2 - cosAngle),2);
				//nowp->lambda /= 1.0/(1.0+Angle);
				//nowp->lambda /= 1.0/(1.0+Angle);
				
				//nowp->lambda /= 1.0/(1.0+Angle);

				//sigsum[i] += sigma[i];
				sigsum[i] += nowp->lambda;
				j++;
			}
			delete [] faceNormal;
			mybcg->sa[i+1] = 1.0;
			mybcg->sa[i+1+numberV] = 1.0;
			vecb[i+1] = 0;			
			vecb[i+1+numberV] = 0;	
		}
		else
		{
			
			mybcg->sa[i+1] = 1.0;
			vecb[i+1] = pU[i];
			mybcg->sa[i+1+numberV] = 1.0;
			vecb[i+1+numberV] = pV[i];		
		}
	}
	mybcg->ija[1] = 2*(numberV)+2;
	int dlk=2*(numberV)+1;
    
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			nowp = PHead[i];
	
			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);
				++dlk;
				mybcg->ija[dlk] = nowp->ID+1;
				mybcg->sa[dlk] = -nowp->lambda/sigsum[i];
				//mybcg->sa[dlk] = -nowp->lambda;

			}
		}
		mybcg->ija[i+1+1]=dlk+1;
	}

	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			nowp = PHead[i];
			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);
				++dlk;
				mybcg->ija[dlk]=nowp->ID+numberV+1;
				mybcg->sa[dlk] = -nowp->lambda/sigsum[i];
				//mybcg->sa[dlk] = -nowp->lambda;
			}
		}
		mybcg->ija[i+numberV+1+1]=dlk+1;
	}
     
	for(i=0;i<numberV;i++)
	{
		UaXY[i+1] = pU[i];
		UaXY[i+numberV+1] = pV[i];
	}

	mybcg->linbcg(((unsigned long)(2*(numberV))),vecb,UaXY,1,PCBCGerror,iteNum,&iter,&linerr);
	for (int i=0;i<numberV;i++)
	{
		pIPV[i].u = pU[i] = UaXY[i+1];
		pIPV[i].v = pV[i] = UaXY[i+numberV+1];		
	}
	
	/*
	SUPERLU_Solve();
	for (int i=0;i<numberV;i++)
	{
		pIPV[i].u = pU[i];
		pIPV[i].v = pV[i];
	}
	*/

	delete [] vecb;
	delete [] UaXY;
	delete [] sigsum; 
	return getCurrentE();
}



void	MyParameterization::setDAMC()
{
	int i;	
	IDList *now=NULL;
	PolarList *nowp=NULL;
	PolarList *nowp_back=NULL;	
	PolarList *nowp_next=NULL;	
	double pi = acos(-1.0);
	PointTool pt;
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
						
			nowp = PHead[i];
			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);	
				//find angle  A , B , C , D
				Point3d *p0 = point[i];
				Point3d *p1 = point[nowp->ID];
				
				nowp_back = nowp->back;
				if (nowp_back ==  PHead[i])
					nowp_back = PTail[i]->back;
				
				nowp_next = nowp->next;
				if (nowp_next ==  PTail[i])
					nowp_next = PHead[i]->next;

				Point3d *p2 = point[nowp_next->ID];
				Point3d *p3 = point[nowp_back->ID];				

				Point3d v0,v1,v2,v3,v4;
				pt.makeVector(&v0,p1,p0);
				pt.makeVector(&v1,p1,p2);
				pt.makeVector(&v2,p1,p3);
				pt.makeVector(&v3,p2,p0);
				pt.makeVector(&v4,p3,p0);
				double Distance = pt.InnerProduct(&v0,&v0);
				pt.Normalize3D(&v0);
				pt.Normalize3D(&v1);
				pt.Normalize3D(&v2);
				pt.Normalize3D(&v3);
				pt.Normalize3D(&v4);
				double AngleA = acos(pt.InnerProduct(&v0,&v1));
				double AngleB = acos(pt.InnerProduct(&v0,&v2));
				pt.ScalarVector(&v1,-1,&v1);
				pt.ScalarVector(&v2,-1,&v2);
				double AngleC = acos(pt.InnerProduct(&v3,&v1));
				double AngleD = acos(pt.InnerProduct(&v4,&v2));
			
				double MX = ((1/tan(AngleA))+(1/tan(AngleB)))/Distance;
				double MA = (1/tan(AngleC))+(1/tan(AngleD));
				nowp->lambda = 0*MA + 1*MX;
				if (nowp->lambda<0)
				{
					int aaa = 0;
				}
			}
		}
	}
}
void	MyParameterization::setPolarMap_EX()
{
	
	int i,j,k;

	VList *now=NULL;
	VList *now2=NULL;
	double dx,dy;
	PolarList *nowp=NULL;
	PolarList *nowp2=NULL;
	int nowID,nextID;
	double angle = 0.0;
	double theta = 0.0;
	int cn0=0;
	// make Geodesic Polar Map 
	for(i=0;i<numberV;i++)
	{
		IDtool->CleanNeighborPolar(PHead[i],PTail[i]);
      
		PHead[i] = new PolarList();
		PTail[i] = new PolarList();
		PHead[i]->next = PTail[i];
		PTail[i]->back = PHead[i];

		if(boundary[i]==0)
		{      
			dy = 0.0;
			theta = 0.0;
			double d =	-(VHead[i]->normal.x*point[i]->x)
						-(VHead[i]->normal.y*point[i]->y)
						-(VHead[i]->normal.z*point[i]->z);
			Point3d *newPT = new Point3d[neighborI[i]];
			now = VHead[i];
			VList *prevV = back(VTail[i]);
			j = 0;
			while(next(now)!=VTail[i])
			{
				now = next(now);
				Point3d N;
				N.x = now->normal.x + prevV->normal.x;
				N.y = now->normal.y + prevV->normal.y;
				N.z = now->normal.z + prevV->normal.z;
				PT->Normalize3D(&N);
				{
					double cosAngle = PT->InnerProduct(&VHead[i]->normal,&N);
					double D = ((VHead[i]->normal.x*point[now->ID]->x)+
								(VHead[i]->normal.y*point[now->ID]->y)+
								(VHead[i]->normal.z*point[now->ID]->z)+d);				
					double L = fabs(D/cosAngle);
					double cos2Angle = 2.0*cosAngle*cosAngle - 1;
					double tanA = sqrt((1.0 - cos2Angle)/(1.0 + cos2Angle));
					Point3d v; 
					PT->makeVector(&v,point[i],point[now->ID]);
					PT->Normalize3D(&v);
					newPT[j].x = point[now->ID]->x + (L*tanA*v.x);
					newPT[j].y = point[now->ID]->y + (L*tanA*v.y);
					newPT[j].z = point[now->ID]->z + (L*tanA*v.z);					
					j = j++;
				}
				now = next(now);
				prevV = now;
			}			
			
			for (j = 0; j < neighborI[i]; j++)
			{
				PT->makeVector(bc[0],point[i],&newPT[j%neighborI[i]]);
				PT->makeVector(bc[1],point[i],&newPT[(j+1)%neighborI[i]]);
				angle = acos((PT->InnerProduct(bc[0],bc[1])/(PT->Point3dSize(bc[0])*PT->Point3dSize(bc[1]))));
				theta+= angle;
			}

			now = VHead[i];
			now = next(now);
			dx = PT->Distance(point[i],&newPT[0]);
			dy = 0.0;
			IDtool->AppendPolarI(now->ID,PTail[i],dx,dy);
			j = 1;
			while(j < neighborI[i])
			{
				now = next(next(now));
				PT->makeVector(bc[1],point[i],&newPT[j%neighborI[i]]);
				PT->makeVector(bc[0],point[i],&newPT[(j-1)%neighborI[i]]);
				dx = PT->Point3dSize(bc[1]);
				angle = acos((PT->InnerProduct(bc[0],bc[1])/(PT->Point3dSize(bc[0])*PT->Point3dSize(bc[1]))));
				dy += (2.0*PI*angle)/theta;
				IDtool->AppendPolarI(now->ID,PTail[i],dx,dy);
				j++;
			}
		}
	}
}

void	MyParameterization::ResetInnerLambda()
{
	for(int i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			PolarList *nowp = PHead[i];
			while(nextN(nowp)!=PTail[i])
			{
				nowp = nextN(nowp);
				nowp->lambda  = 0;
				/*
				nowp->lambdaU = 0;
				nowp->lambdaV = 0;
				*/
				
			}
		}
	}
}

double    MyParameterization::PARAM_MYEXPER2(PolarVertex *pIPV,
											 int num_PV,
											 FILE* logFile)
{
	double bestStretch = -1;
	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;	
	if (pU)
		delete [] pU;
	if (pV)
		delete [] pV;

	pU = new double[numberV];
	pV = new double[numberV];
	memset(pU,0x00,sizeof(double)*numberV);
	memset(pV,0x00,sizeof(double)*numberV);
	


    boundarytype=0;
	setPolarMap();
	IDList *BpointH = new IDList();
	IDList *BpointT = new IDList();
    
	BpointH->next = BpointT;
	BpointT->back = BpointH;
	
	double tlength=0.0;
	int numBorderPoint = 0;
	CalBorderPath(BpointH,BpointT,&tlength,&numBorderPoint);

	IDList *now = BpointH;
	
	int state=0;
	double clen=0.0;		
	double loop = 0;
	double sum_length = 0;
	double edge_length = 0;		
	IDList *startPoint = BpointH->next;
		
	int count = 1;
	while (loop < tlength*0.25)
	{
		BpointT->ID = BpointH->next->ID; // for cal length;
		state = 0;
		sum_length = 0;
		double this_side_length[4] ={0};
		int this_side_num_edge[4] = {0};		
		loop += PT->Distance(point[startPoint->ID],point[next(startPoint)->ID]);
		IDList *now = startPoint;
		//if (count == 1)
		{
		do
		{				
			edge_length = PT->Distance(point[now->ID],point[next(now)->ID]);			
			if ((sum_length + edge_length>= (tlength*0.25)) && state < 3)
			{
				if ((sum_length + edge_length)/(tlength*0.25) >= 1.1&& this_side_num_edge[state] > 0)					
				{
					this_side_length[state] = sum_length;
					this_side_num_edge[state+1]++;
					sum_length = edge_length;
				}
				else
				{
					this_side_length[state] = sum_length+edge_length;
					this_side_num_edge[state]++;
					sum_length = 0;
				}				
				state++;
			}
			else
			{
				this_side_num_edge[state]++;
				sum_length += edge_length;
			}
			now = next(now);
			if (now == BpointT)
				now = BpointH->next;
		}
		while (now != startPoint);			
		this_side_length[state] = sum_length;


		state = 0;
		sum_length  = 0;
		clen = 0.0;
		now = BpointH;
		BpointT->ID = BpointH->next->ID; // for cal length;
		while(next(now)!=BpointT)
		{
			now = next(now);
			switch(state)
			{
				case 0:
					if (clen < 1.0)
					{
						pU[now->ID] = clen; 
						pV[now->ID] = 0.0;
					}
					else
					{
						pU[now->ID] = 1.0;
						pV[now->ID] = 0.0;
						sum_length = 0.0;
						state=1;
						
					}
					break;
				case 1:
					if (clen < 1.0)
					{
						pU[now->ID] = 1.0; 
						pV[now->ID] = clen;
					}
					else
					{
						pU[now->ID] = 1.0;
						pV[now->ID] = 1.0;
						sum_length = 0.0;
						state=2;

					}
					break;
				case 2:
					if (clen < 1.0)
					{
						pU[now->ID] = 1.0-clen; 
						pV[now->ID] = 1.0;
					}
					else
					{
						pU[now->ID] = 0.0;
						pV[now->ID] = 1.0;
						sum_length = 0.0;
						state=3;

					}
					break;
				case 3:
					if (clen < 1.0)
					{
						pU[now->ID] = 0.0; 
						pV[now->ID] = 1.0-clen;
					}
					else
					{
						//should never enter this one;
						state=4;
					}
				default:
					break;

			}		
			sum_length += PT->Distance(point[now->ID],point[next(now)->ID]);
			clen = sum_length/this_side_length[state];	
		}
		
		for(int i=0;i<numberV;i++)
		{
			if(boundary[i]==1)
			{            
		
			}
			else
			{
				pU[i] = 0;
				pV[i] = 0;
			}
		}
		ResetInnerLambda();	
		printf("TEST %d\n",count);
		fprintf(logFile,"TEST %d\n",count);
		gammaP = 0.5;
		ParametrizationOptimal(iteNum,PCBCGerror,logFile);
		//ParametrizationSmoothOptimal_EX(iteNum,PCBCGerror,logFile);
		//ParametrizationSmoothOptimal(iteNum,PCBCGerror,logFile);
		if (resultStretch < bestStretch || bestStretch < 0)
		{
			for (int i=0;i<numberV;i++)
			{
				pIPV[i].u = pU[i];
				pIPV[i].v = pV[i];		
			}
			bestStretch = resultStretch;
		}
		}

		if (loop >= tlength*0.25)
			break;
		count++;
		//prepare for next loop
		startPoint = next(startPoint);
		//reconstruct BpointH , BpointT according to startPoint
		int startPointID = startPoint->ID;
		while (BpointH->next->ID != startPoint->ID)
		{
			IDList *move_node = BpointH->next;
			BpointH->next = move_node->next;
			move_node->next->back = BpointH; 
			
			move_node->back = BpointT->back;
			BpointT->back->next = move_node;			
			BpointT->back = move_node;
			move_node->next = BpointT;			
		}
	}
	IDtool->CleanNeighbor(BpointH,BpointT);
	return bestStretch;
}

double    MyParameterization::PARAM_MYEXPER3(PolarVertex *pIPV,
											int num_PV,
											FILE* logFile)
{
	double bestStretch = -1;
	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;	
	if (pU)
		delete [] pU;
	if (pV)
		delete [] pV;

	pU = new double[numberV];
	pV = new double[numberV];
	memset(pU,0x00,sizeof(double)*numberV);
	memset(pV,0x00,sizeof(double)*numberV);
	
	double *initU0 = new double[numberV];
	double *initV0 = new double[numberV];
	memset(initU0,0x00,sizeof(double)*numberV);
	memset(initV0,0x00,sizeof(double)*numberV);

	clock_t start,finish;
	start = clock();

    boundarytype=0;
	setPolarMap();
	IDList *BpointH = new IDList();
	IDList *BpointT = new IDList();
    
	BpointH->next = BpointT;
	BpointT->back = BpointH;
	
	double tlength=0.0;
	int numBorderPoint = 0;
	CalBorderPath(BpointH,BpointT,&tlength,&numBorderPoint);

	IDList *now = BpointH;
	
	int state=0;
	double clen=0.0;		
	double loop = 0;
	double sum_length = 0;
	double edge_length = 0;		
	IDList *startPoint = BpointH->next;
	int best_startID = startPoint->ID;	
	int worst_startID = startPoint->ID;
	double worstStretch = -1;
	int count = 0;

	FILE *fpCorner = fopen("cornerID_setting.txt","rt");
	int spec_cornerID = -1;

	if (fpCorner)
	{
		printf("FOUND CORNER SETTING FILE!\n");
		//read setting from file		
		fscanf(fpCorner,"%d\r\n",&spec_cornerID);		
		fclose(fpCorner);
		if (spec_cornerID >= 0)
			printf("Calculate CORNER ID%d\n",spec_cornerID);
		else
		{
			spec_cornerID = -1;
			printf("Calculate ALL corner\n");
		}
		 
	}

	FILE *fpCornerTestResult = fopen("corner_length_stretch.txt","wt");
	while (loop < tlength*0.25)
	{
		BpointT->ID = BpointH->next->ID; // for cal length;
		state = 0;
		sum_length = 0;
		double this_side_length[4] ={0};
		int this_side_num_edge[4] = {0};
		fprintf(fpCornerTestResult,"%f\t",loop);
		loop += PT->Distance(point[startPoint->ID],point[next(startPoint)->ID]);
		IDList *now = startPoint;
		if (spec_cornerID < 0 || spec_cornerID == count)
		{
			do
			{				
				edge_length = PT->Distance(point[now->ID],point[next(now)->ID]);			
				if ((sum_length + edge_length>= (tlength*0.25)) && state < 3)
				{
					if ((sum_length + edge_length)/(tlength*0.25) >= 1.1&& this_side_num_edge[state] > 0)					
					{
						this_side_length[state] = sum_length;
						this_side_num_edge[state+1]++;
						sum_length = edge_length;
					}
					else
					{
						this_side_length[state] = sum_length+edge_length;
						this_side_num_edge[state]++;
						sum_length = 0;
					}				
					state++;
				}
				else
				{
					this_side_num_edge[state]++;
					sum_length += edge_length;
				}
				now = next(now);
				if (now == BpointT)
					now = BpointH->next;
			}
			while (now != startPoint);			
			this_side_length[state] = sum_length;

			//memcpy(pU,initU0,sizeof(double)*numberV);
			//memcpy(pV,initV0,sizeof(double)*numberV);

			state = 0;
			sum_length  = 0;
			clen = 0.0;
			now = BpointH;
			BpointT->ID = BpointH->next->ID; // for cal length;
		
			printf("TEST %d (%f/%f)  POINT CORNER ID: ",count,loop,tlength*0.25);		
			printf("%d ",BpointH->next->ID);
			
			if (logFile)
			{
				fprintf(logFile,"TEST %d (%f/%f)  POINT CORNER ID: ",count,loop,tlength*0.25);
				fprintf(logFile,"%d ",BpointH->next->ID);
			}
			
				
			while(next(now)!=BpointT)
			{
				now = next(now);
				switch(state)
				{
					case 0:
						if (clen < 1.0)
						{
							pU[now->ID] = clen; 
							pV[now->ID] = 0.0;
						}
						else
						{
							pU[now->ID] = 1.0;
							pV[now->ID] = 0.0;
							sum_length = 0.0;
							state=1;
							printf("%d ",now->ID);
							
							if (logFile)
								fprintf(logFile,"%d ",now->ID);
							
						}
						break;
					case 1:
						if (clen < 1.0)
						{
							pU[now->ID] = 1.0; 
							pV[now->ID] = clen;
						}
						else
						{
							pU[now->ID] = 1.0;
							pV[now->ID] = 1.0;
							sum_length = 0.0;
							state=2;
							printf("%d ",now->ID);
							
							if (logFile)
								fprintf(logFile,"%d ",now->ID);
							

						}
						break;
					case 2:
						if (clen < 1.0)
						{
							pU[now->ID] = 1.0-clen; 
							pV[now->ID] = 1.0;
						}
						else
						{
							pU[now->ID] = 0.0;
							pV[now->ID] = 1.0;
							sum_length = 0.0;
							state=3;
							printf("%d\n",now->ID);
							
							if (logFile)
								fprintf(logFile,"%d\n",now->ID);
							

						}
						break;
					case 3:
						if (clen < 1.0)
						{
							pU[now->ID] = 0.0; 
							pV[now->ID] = 1.0-clen;
						}
						else
						{
							//should never enter this one;
							state=4;
						}
					default:
						break;

				}		
				sum_length += PT->Distance(point[now->ID],point[next(now)->ID]);
				clen = sum_length/this_side_length[state];	
			}
		
			for(int i=0;i<numberV;i++)
			{
				if(boundary[i]==1)
				{            
		
				}
				else
				{
					pU[i] = 0.0;
					pV[i] = 0.0;
				}
			}
		
			ResetInnerLambda();		
			gammaP = 0.75;
			//ParametrizationOptimalSaveU0(initU0,initV0,iteNum,PCBCGerror,NULL);
			ParametrizationOptimal(iteNum,PCBCGerror,NULL);
			fprintf(fpCornerTestResult,"%f\n",resultStretch);
			if (logFile)
			{
				fprintf(logFile,"STRETCH=%f\n",resultStretch);
				
			}
			gammaP = 1.0;
			if (resultStretch < bestStretch || bestStretch < 0)
			{
				for (int i=0;i<numberV;i++)
				{
					if(boundary[i]!=1)
					{
						PolarList *nowp = PHead[i];
						while(next(nowp)!=PTail[i])
						{
							nowp = next(nowp);
							nowp->best_lambda = nowp->old_lambda;
						}
					}
					pIPV[i].u = pU[i];
					pIPV[i].v = pV[i];
				}
				bestStretch = resultStretch;
				best_startID = count;
			}

			if (resultStretch > worstStretch)
			{
				worstStretch = resultStretch;
				worst_startID = count;
			}
		}

		if (loop >= tlength*0.25)
			break;

		count++;
		//prepare for next loop
		startPoint = next(startPoint);
		//reconstruct BpointH , BpointT according to startPoint
		int startPointID = startPoint->ID;
		while (BpointH->next->ID != startPoint->ID)
		{
			IDList *move_node = BpointH->next;
			BpointH->next = move_node->next;
			move_node->next->back = BpointH;
			
			move_node->back = BpointT->back;
			BpointT->back->next = move_node;
			BpointT->back = move_node;
			move_node->next = BpointT;
		}
	}
	fclose(fpCornerTestResult);
	delete [] initU0;
	delete [] initV0;


	//reconstruct BpointH , BpointT according to startPoint	
	/*
	while (BpointH->next->ID != best_startID)
	{
		IDList *move_node = BpointH->next;
		BpointH->next = move_node->next;
		move_node->next->back = BpointH;
			
		move_node->back = BpointT->back;
		BpointT->back->next = move_node;
		BpointT->back = move_node;
		move_node->next = BpointT;
	}
	*/


	finish = clock();
	printf( "Calculate Time for best corner : %f sec\n", (double)(finish-start)/CLOCKS_PER_SEC);
	if (logFile)
		fprintf(logFile,"Calculate Time for best corner : %f sec\n", (double)(finish-start)/CLOCKS_PER_SEC);
	

	for (int i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			PolarList *nowp = PHead[i];
			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);
				nowp->lambda = nowp->best_lambda;
			}
		}
		pU[i] = pIPV[i].u;
		pV[i] = pIPV[i].v;
	}

	printf("=== Find best stretch  of TEST%d  (best corner at %f) ===\n",best_startID,bestStretch);
	if (logFile)
		fprintf(logFile,"=== Find best stretch of TEST%d (best corner at %f)===\n",best_startID,bestStretch);

	ParametrizationSmoothOptimal_EX(0.5,0.1,iteNum,PCBCGerror,logFile);

	printf("=== Find worst stretch  of TEST%d  (worst corner at %f) ===\n",worst_startID,worstStretch);
	if (logFile)
		fprintf(logFile,"=== Find worst stretch  of TEST%d  (worst corner at %f) ===\n",worst_startID,worstStretch);
	
	for (int i=0;i<numberV;i++)
	{
		pIPV[i].u = pU[i];
		pIPV[i].v = pV[i];
	}
	IDtool->CleanNeighbor(BpointH,BpointT);
	return resultStretch;
}


//search best corner point by local minimum search
double    MyParameterization::PARAM_MYEXPER4(PolarVertex *pIPV,
											int num_PV,
											FILE* logFile)
{
	double bestStretch = -1;
	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;	
	if (pU)
		delete [] pU;
	if (pV)
		delete [] pV;

	pU = new double[numberV];
	pV = new double[numberV];
	memset(pU,0x00,sizeof(double)*numberV);
	memset(pV,0x00,sizeof(double)*numberV);
	
	double *initU0 = new double[numberV];
	double *initV0 = new double[numberV];
	memset(initU0,0x00,sizeof(double)*numberV);
	memset(initV0,0x00,sizeof(double)*numberV);

	clock_t start,finish;
	start = clock();

    boundarytype=0;
	setPolarMap();
	IDList *BpointH = new IDList();
	IDList *BpointT = new IDList();
    
	BpointH->next = BpointT;
	BpointT->back = BpointH;
	
	double tlength=0.0;
	int numBorderPoint = 0;
	CalBorderPath(BpointH,BpointT,&tlength,&numBorderPoint);

	IDList *now = BpointH;
	
	int state=0;
	double clen=0.0;		
	double loop = 0;
	double sum_length = 0;
	double edge_length = 0;		
	IDList *startPoint = BpointH->next;
	int best_startID = startPoint->ID;	
	int worst_startID = startPoint->ID;
	double worstStretch = -1;
	int count = 0;

	FILE *fpCorner = fopen("cornerID_setting.txt","rt");
	int spec_cornerID = -1;

	if (fpCorner)
	{
		printf("FOUND CORNER SETTING FILE!\n");
		//read setting from file		
		fscanf(fpCorner,"%d\r\n",&spec_cornerID);		
		fclose(fpCorner);
		if (spec_cornerID >= 0)
			printf("Calculate CORNER ID%d\n",spec_cornerID);
		else
		{
			spec_cornerID = -1;
			printf("Calculate ALL corner\n");
		}
		 
	}

	//record result as array.
	double *stretch_each_corner = new double[numBorderPoint];
	while (loop < tlength*0.25)
	{
		BpointT->ID = BpointH->next->ID; // for cal length;
		state = 0;
		sum_length = 0;
		double this_side_length[4] ={0};
		int this_side_num_edge[4] = {0};		
		loop += PT->Distance(point[startPoint->ID],point[next(startPoint)->ID]);
		IDList *now = startPoint;
		if (spec_cornerID < 0 || spec_cornerID == count)
		{
			do
			{				
				edge_length = PT->Distance(point[now->ID],point[next(now)->ID]);			
				if ((sum_length + edge_length>= (tlength*0.25)) && state < 3)
				{
					if ((sum_length + edge_length)/(tlength*0.25) >= 1.1&& this_side_num_edge[state] > 0)					
					{
						this_side_length[state] = sum_length;
						this_side_num_edge[state+1]++;
						sum_length = edge_length;
					}
					else
					{
						this_side_length[state] = sum_length+edge_length;
						this_side_num_edge[state]++;
						sum_length = 0;
					}				
					state++;
				}
				else
				{
					this_side_num_edge[state]++;
					sum_length += edge_length;
				}
				now = next(now);
				if (now == BpointT)
					now = BpointH->next;
			}
			while (now != startPoint);			
			this_side_length[state] = sum_length;

			//memcpy(pU,initU0,sizeof(double)*numberV);
			//memcpy(pV,initV0,sizeof(double)*numberV);

			state = 0;
			sum_length  = 0;
			clen = 0.0;
			now = BpointH;
			BpointT->ID = BpointH->next->ID; // for cal length;
		
			printf("TEST %d (%f/%f)  POINT CORNER ID: ",count,loop,tlength*0.25);		
			printf("%d ",BpointH->next->ID);
			
			if (logFile)
			{
				fprintf(logFile,"TEST %d (%f/%f)  POINT CORNER ID: ",count,loop,tlength*0.25);
				fprintf(logFile,"%d ",BpointH->next->ID);
			}
			
				
			while(next(now)!=BpointT)
			{
				now = next(now);
				switch(state)
				{
					case 0:
						if (clen < 1.0)
						{
							pU[now->ID] = clen; 
							pV[now->ID] = 0.0;
						}
						else
						{
							pU[now->ID] = 1.0;
							pV[now->ID] = 0.0;
							sum_length = 0.0;
							state=1;
							printf("%d ",now->ID);
							
							if (logFile)
								fprintf(logFile,"%d ",now->ID);
							
						}
						break;
					case 1:
						if (clen < 1.0)
						{
							pU[now->ID] = 1.0; 
							pV[now->ID] = clen;
						}
						else
						{
							pU[now->ID] = 1.0;
							pV[now->ID] = 1.0;
							sum_length = 0.0;
							state=2;
							printf("%d ",now->ID);
							
							if (logFile)
								fprintf(logFile,"%d ",now->ID);
							

						}
						break;
					case 2:
						if (clen < 1.0)
						{
							pU[now->ID] = 1.0-clen; 
							pV[now->ID] = 1.0;
						}
						else
						{
							pU[now->ID] = 0.0;
							pV[now->ID] = 1.0;
							sum_length = 0.0;
							state=3;
							printf("%d\n",now->ID);
							
							if (logFile)
								fprintf(logFile,"%d\n",now->ID);
							

						}
						break;
					case 3:
						if (clen < 1.0)
						{
							pU[now->ID] = 0.0; 
							pV[now->ID] = 1.0-clen;
						}
						else
						{
							//should never enter this one;
							state=4;
						}
					default:
						break;

				}		
				sum_length += PT->Distance(point[now->ID],point[next(now)->ID]);
				clen = sum_length/this_side_length[state];	
			}
		
			for(int i=0;i<numberV;i++)
			{
				if(boundary[i]==1)
				{            
		
				}
				else
				{
					pU[i] = 0.0;
					pV[i] = 0.0;
				}
			}

			ResetInnerLambda();		
			gammaP = 0.75;
			//ParametrizationOptimalSaveU0(initU0,initV0,iteNum,PCBCGerror,NULL);
			ParametrizationOptimal(iteNum,PCBCGerror,NULL);
			printf("STRETCH=%f\n",resultStretch);
			if (logFile)
			{
				fprintf(logFile,"STRETCH=%f\n",resultStretch);
				
			}
			gammaP = 1.0;
			if (resultStretch < bestStretch || bestStretch < 0)
			{
				for (int i=0;i<numberV;i++)
				{
					if(boundary[i]!=1)
					{
						PolarList *nowp = PHead[i];
						while(next(nowp)!=PTail[i])
						{
							nowp = next(nowp);
							nowp->best_lambda = nowp->old_lambda;
						}
					}
					pIPV[i].u = pU[i];
					pIPV[i].v = pV[i];
				}
				bestStretch = resultStretch;
				best_startID = count;
			}

			if (resultStretch > worstStretch)
			{
				worstStretch = resultStretch;
				worst_startID = count;
			}
		}

		if (loop >= tlength*0.25)
			break;

		count++;
		//prepare for next loop
		startPoint = next(startPoint);
		//reconstruct BpointH , BpointT according to startPoint
		int startPointID = startPoint->ID;
		while (BpointH->next->ID != startPoint->ID)
		{
			IDList *move_node = BpointH->next;
			BpointH->next = move_node->next;
			move_node->next->back = BpointH;
			
			move_node->back = BpointT->back;
			BpointT->back->next = move_node;
			BpointT->back = move_node;
			move_node->next = BpointT;
		}
	}
	
	delete [] initU0;
	delete [] initV0;


	finish = clock();
	printf( "Calculate Time for best corner : %f sec\n", (double)(finish-start)/CLOCKS_PER_SEC);
	if (logFile)
		fprintf(logFile,"Calculate Time for best corner : %f sec\n", (double)(finish-start)/CLOCKS_PER_SEC);
	
	/*
	for (int i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			PolarList *nowp = PHead[i];
			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);
				nowp->lambda = nowp->best_lambda;
			}
		}
		pU[i] = pIPV[i].u;
		pV[i] = pIPV[i].v;
	}

	printf("=== Find best stretch  of TEST%d  (best corner at %f) ===\n",best_startID,bestStretch);
	if (logFile)
		fprintf(logFile,"=== Find best stretch of TEST%d (best corner at %f)===\n",best_startID,bestStretch);

	ParametrizationSmoothOptimal_EX(0.5,0.1,iteNum,PCBCGerror,logFile);

	printf("=== Find worst stretch  of TEST%d  (worst corner at %f) ===\n",worst_startID,worstStretch);
	if (logFile)
		fprintf(logFile,"=== Find worst stretch  of TEST%d  (worst corner at %f) ===\n",worst_startID,worstStretch);
	
	for (int i=0;i<numberV;i++)
	{
		pIPV[i].u = pU[i];
		pIPV[i].v = pV[i];
	}
	*/
	IDtool->CleanNeighbor(BpointH,BpointT);
	return resultStretch;
}


double    MyParameterization::PARAM_MYEXPER(PolarVertex *pIPV,
										int num_PV,
										FILE* logFile)
{
	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;	
    boundarytype=0;
	MyBoundaryMap();
	//BoundaryMap();
	setPolarMap();


	//setPolarMap_EX();

	double *UaXY = new double[2*(numberV)+1];
    double *vecb = new double[2*(numberV)+1];
  
	int i;
	IDList *now;
	IDList *now2;
	PolarList *nowp;
	level=0;
  
	int nonzero=(numberV);
	for(i=0;i<numberV;i++)
	{
		vecb[i+1]=0.0;
		if(boundary[i]!=1)
		{
			now = IHead[i];
			while(next(now)!=ITail[i])
			{
				now = next(now);
				nonzero++;
			}
		}
	}

	int iter=0;
	double linerr=0.0;
	double weight=0.0;
  
	PCBCGSolver *mybcg = new PCBCGSolver(2*nonzero);
	double *sigsum = new double[numberV];
	double *sigsumU = new double[numberV];
	double *sigsumV = new double[numberV];
	setFloaterC();

    
  
	SortIndexP();
  
  

	for(i=0;i<numberV;i++)
	{    
		if(boundary[i]!=1)
		{
			mybcg->sa[i+1] = 1.0;
			mybcg->sa[i+1+numberV] = 1.0;
			vecb[i+1] = 0.0;
			vecb[i+1+numberV] = 0.0;
		}
		else
		{
			mybcg->sa[i+1] = 1.0;
			vecb[i+1] = pU[i];
			mybcg->sa[i+1+numberV] = 1.0;
			vecb[i+1+numberV] = pV[i];
		}
	}

	mybcg->ija[1] = 2*(numberV)+2;
	int dlk=2*(numberV)+1;
  
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			nowp = PHead[i];      
			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);
				++dlk;
				nowp->lambdaU = nowp->lambda;
				nowp->lambdaV = nowp->lambda;
				mybcg->sa[dlk] = -nowp->lambda;
				mybcg->ija[dlk]=nowp->ID+1;
			}
		}
		mybcg->ija[i+1+1]=dlk+1;
	}
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			nowp = PHead[i];
			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);
				++dlk;
				mybcg->sa[dlk] = -nowp->lambda;
				mybcg->ija[dlk]=nowp->ID+numberV+1;
			}
		}
		mybcg->ija[i+numberV+1+1]=dlk+1;
	}
  
	for(i=0;i<numberV;i++)
	{
		UaXY[i+1] = pU[i];
		UaXY[i+numberV+1] = pV[i];
	}
	mybcg->linbcg(((unsigned long)(2*(numberV))),vecb,UaXY,1,PCBCGerror,iteNum,&iter,&linerr);
  
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			pU[i] = UaXY[i+1];
			pV[i] = UaXY[i+numberV+1];
		}
	}
  
	// Re-solving linear system
	//double initialstrech=0.0;
	//double currentstrech=0.0;
	double initialstrechU=0.0;
	double currentstrechU=0.0;
	double initialstrechV=0.0;
	double currentstrechV=0.0;
	Point2d *prevUV = new Point2d[numberV];	

	initialstrechU = getCurrentE_U();
	int kk = 0;  

	printf("R%d  STRETCH U: %f\n",kk,initialstrechU); //u0
	if (logFile)
	{
		fprintf(logFile,"R%d  STRETCH U: %f\n",kk,initialstrechU); //u0
	}
	initialstrechV = getCurrentE_V();
	printf("R%d  STRETCH V: %f\n",kk,initialstrechV); //u0
	if (logFile)
	{
		fprintf(logFile,"R%d  STRETCH V: %f\n",kk,initialstrechV); //u0
	}

	//solve for UV
	//for(kk=0;kk<iteNum;kk++)
	{    
		
		while (gammaU >= 1.0/16.0)	
		{
			setSigmaU();
			//U section
			for(i=0;i<numberV;i++)
			{      
				if(boundary[i]!=1)
				{
					//sigsum[i]=0.0;
					sigsumU[i] = 0.0;					
					nowp = PHead[i];
					while(next(nowp)!=PTail[i])
					{
						nowp = next(nowp);
						nowp->old_lambdaU = nowp->lambdaU;						
						nowp->lambdaU /= sigmaU[nowp->ID];						
						sigsumU[i] += nowp->lambdaU;						
					}
					mybcg->sa[i+1] = 1.0;
					
				}
				else
				{
					mybcg->sa[i+1] = 1.0;
					
				}
			}
			dlk=(numberV)+1;
    
			for(i=0;i<numberV;i++)
			{
				if(boundary[i]!=1)
				{
					nowp = PHead[i];
	
					while(next(nowp)!=PTail[i])
					{
						nowp = next(nowp);
						++dlk;
						//mybcg->sa[dlk] = -nowp->lambda/sigsum[i];
						mybcg->sa[dlk] = -nowp->lambdaU/sigsumU[i];
						mybcg->ija[dlk] = nowp->ID+1;
					}
				}
				mybcg->ija[i+1+1]=dlk+1;
			}
		
			for(i=0;i<numberV;i++)
			{
				UaXY[i+1] = pU[i];
				//UaXY[i+numberV+1] = pV[i];
			}

			//mybcg->linbcg(((unsigned long)(2*(numberV))),vecb,UaXY,1,PCBCGerror,iteNum,&iter,&linerr); 
			mybcg->ija[1] = (numberV)+2;
			mybcg->linbcg(((unsigned long)((numberV))),vecb,UaXY,1,PCBCGerror,iteNum,&iter,&linerr);
			bool failToSolve = false;
			if (iter > iteNum)
			{
				failToSolve = true;				
			}
		
			for(i=0;i<numberV;i++)
			{
				if(boundary[i]!=1)
				{
					prevUV[i].x = pU[i];
					//prevU[i]->y = pV[i];
	
					pU[i] = UaXY[i+1];
					//pV[i] = UaXY[i+numberV+1];
				}
			}
			currentstrechU = getCurrentE_U();
			//printf("currentstrech U%d= %lf\n",kk+1,currentstrech);
			printf("R%d  STRETCH U: %f\n",kk+1,currentstrechU); 

			if (logFile)
			{
				fprintf(logFile,"R%d  STRETCH U: %f\n",kk+1,currentstrechU); 
			}

			if(initialstrechU<=currentstrechU || failToSolve)
			{
			
				for(i=0;i<numberV;i++)
				{
					if(boundary[i]!=1)
					{
						pU[i] = prevUV[i].x;					
						nowp = PHead[i];
						while(next(nowp)!=PTail[i])
						{
							nowp = next(nowp);
							nowp->lambdaU = nowp->old_lambdaU;
						}
					}
				}
				gammaU /= 2.0;
			}
			else
			{
				initialstrechU = currentstrechU;
			}
		}
		
		
		while (gammaV >= 1.0/16.0)	
		{
			setSigmaV();
			//V section
			for(i=0;i<numberV;i++)
			{      
				if(boundary[i]!=1)
				{
					//sigsum[i]=0.0;
					//sigsumU[i] = 0.0;
					  
  
					sigsumV[i] = 0.0;
					nowp = PHead[i];
					while(next(nowp)!=PTail[i])
					{
						nowp = next(nowp);					
						nowp->old_lambdaV = nowp->lambdaV ; 
						nowp->lambdaV /= sigmaV[nowp->ID];					
						sigsumV[i] += nowp->lambdaV;
					}
					mybcg->sa[i+1] = 1.0;
					
				}
				else
				{
					mybcg->sa[i+1] = 1.0;
					
				}
			}
			dlk=(numberV)+1;
    
			for(i=0;i<numberV;i++)
			{
				if(boundary[i]!=1)
				{
					nowp = PHead[i];
	
					while(next(nowp)!=PTail[i])
					{
						nowp = next(nowp);
						++dlk;						
						mybcg->sa[dlk] = -nowp->lambdaV/sigsumV[i];
						mybcg->ija[dlk] = nowp->ID+1;
					}
				}
				mybcg->ija[i+1+1]=dlk+1;
			}		
			for(i=0;i<numberV;i++)
			{
				UaXY[i+1] = pV[i];				
			}

			//mybcg->linbcg(((unsigned long)(2*(numberV))),vecb,UaXY,1,PCBCGerror,iteNum,&iter,&linerr); 
			mybcg->ija[1] = (numberV)+2;
			mybcg->linbcg(((unsigned long)((numberV))),vecb,UaXY,1,PCBCGerror,iteNum,&iter,&linerr); 
			bool failToSolve = false;
			if (iter > iteNum)
			{
				failToSolve = true;
			}
		
			for(i=0;i<numberV;i++)
			{
				if(boundary[i]!=1)
				{
					prevUV[i].y = pV[i];	
					pV[i] = UaXY[i+1];
				}
			}
			currentstrechV = getCurrentE_V();
			//printf("currentstrech U%d= %lf\n",kk+1,currentstrech);
			printf("R%d  STRETCH V: %f\n",kk+1,currentstrechV); 

			if (logFile)
			{
				fprintf(logFile,"R%d  STRETCH V: %f\n",kk+1,currentstrechV); 
			}

			if(initialstrechV<=currentstrechV || failToSolve)
			{
			
				for(i=0;i<numberV;i++)
				{
					if(boundary[i]!=1)
					{
						pV[i] = prevUV[i].y;
						nowp = PHead[i];
						while(next(nowp)!=PTail[i])
						{
							nowp = next(nowp);
							nowp->lambdaV = nowp->old_lambdaV;
						}
					}
				}
				gammaV /= 2.0;			
			}
			else
			{
				initialstrechV = currentstrechV;
			}
		}
		
		//if (gammaU < 1.0/16.0 && gammaV < 1.0/16.0)
		//	break;
		
	}


	delete mybcg;	
	
	
	delete [] prevUV;	
	
	delete [] vecb;
	delete [] UaXY;
	delete [] sigsum;
	delete [] sigsumU;
	delete [] sigsumV;
	
	for (int i=0;i<numberV;i++)
	{
		pIPV[i].u = pU[i];
		pIPV[i].v = pV[i];
	}	
	return 0;
}

void MyParameterization::createMirrorFaceData()
{
	mirrorFace = new MirrorFace[numberF];
	PointTool pt;
	for (int i = 0 ; i < numberF ; i++)
	{
		/*if (boundary[Face[i][0]] == 0 &&
			boundary[Face[i][1]] == 0 &&
			boundary[Face[i][2]] == 0)*/
		{
			for (int j = 0 ; j < 3 ; j++)
			{
				int pt_id = Face[i][j];
				//pt.Add(&mirrorFace[i].pt[j],point[pt_id]);
				VList *nowv = VHead[pt_id];
				double d =	-(VHead[pt_id]->normal.x*point[pt_id]->x)
							-(VHead[pt_id]->normal.y*point[pt_id]->y)
							-(VHead[pt_id]->normal.z*point[pt_id]->z);
				while (next(nowv) != VTail[pt_id])
				{
					nowv = next(nowv);
					if (nowv->FaceID == i)
					{
						
						double D = ((VHead[pt_id]->normal.x*point[nowv->ID]->x)+
									(VHead[pt_id]->normal.y*point[nowv->ID]->y)+
									(VHead[pt_id]->normal.z*point[nowv->ID]->z)+d);
						double cosAngle = pt.InnerProduct(&VHead[pt_id]->normal,&nowv->normal);
						if (cosAngle > 1)
							cosAngle = 1;
						if (cosAngle < -1)
							cosAngle = -1;
						double L = fabs(D/cosAngle);
						double cos2Angle = 2.0*cosAngle*cosAngle - 1;
						double tanA = sqrt((1.0 - cos2Angle)/(1.0 + cos2Angle));						


						Point3d newPT;
						Point3d v; 
						pt.makeVector(&v,point[pt_id],point[nowv->ID]);
						pt.Normalize3D(&v);

						newPT.x = point[nowv->ID]->x + (L*tanA*v.x);
						newPT.y = point[nowv->ID]->y + (L*tanA*v.y);
						newPT.z = point[nowv->ID]->z + (L*tanA*v.z);

						if (Face[i][(j+1)%3] == nowv->ID)
							pt.Add(&mirrorFace[i].pt[(j+1)%3],&newPT);
						else
							pt.Add(&mirrorFace[i].pt[(j+2)%3],&newPT);
					}
					nowv = next(nowv);
					if (nowv->FaceID == i)
					{						
						double D = ((VHead[pt_id]->normal.x*point[nowv->ID]->x)+
									(VHead[pt_id]->normal.y*point[nowv->ID]->y)+
									(VHead[pt_id]->normal.z*point[nowv->ID]->z)+d);
						double cosAngle = pt.InnerProduct(&VHead[pt_id]->normal,&nowv->normal);
						if (cosAngle > 1)
							cosAngle = 1;
						if (cosAngle < -1)
							cosAngle = -1;
						double L = fabs(D/cosAngle);
						double cos2Angle = 2.0*cosAngle*cosAngle - 1;
						double tanA = sqrt((1.0 - cos2Angle)/(1.0 + cos2Angle));						
						

						Point3d newPT;
						Point3d v; 
						pt.makeVector(&v,point[pt_id],point[nowv->ID]);
						pt.Normalize3D(&v);

						newPT.x = point[nowv->ID]->x + (L*tanA*v.x);
						newPT.y = point[nowv->ID]->y + (L*tanA*v.y);
						newPT.z = point[nowv->ID]->z + (L*tanA*v.z);

						if (Face[i][(j+1)%3] == nowv->ID)
							pt.Add(&mirrorFace[i].pt[(j+1)%3],&newPT);
						else
							pt.Add(&mirrorFace[i].pt[(j+2)%3],&newPT);

						break;
					}
				}
			}
			for (int j = 0 ; j < 3 ; j++)
			{
				pt.ScalarVector(&mirrorFace[i].pt[j],1.0/2.0,&mirrorFace[i].pt[j]);
			}
			mirrorFace[i].area3d = pt.getArea(&mirrorFace[i].pt[0],&mirrorFace[i].pt[1],&mirrorFace[i].pt[2]);			

			if (mirrorFace[i].area3d < areaMap3D[i])
				int aaa = 0;
		}
		/*
		else
		{
			pt.Copy(point[Face[i][0]],  &mirrorFace[i].pt[0]);
			pt.Copy(point[Face[i][1]],  &mirrorFace[i].pt[1]);
			pt.Copy(point[Face[i][2]],  &mirrorFace[i].pt[2]);
			mirrorFace[i].area3d = areaMap3D[i];
		}
		*/
	}
}


void MyParameterization::setSigmaArea()
{
	
	int i,j;
	IDList *now=NULL;
	double varphi,ddv,dsize1,sumarea;
	double dddhval=0.0;
	double localsum=0.0;
	for(i=0;i<numberF;i++)
	{
    
		dsize1 = PT->getParametricA(pV[Face[i][0]],
					pV[Face[i][1]],
					pV[Face[i][2]],
					pU[Face[i][0]],
					pU[Face[i][1]],
					pU[Face[i][2]]);
		
		PT->setParametricDs(bc[0],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pV[Face[i][0]],pV[Face[i][1]],
				pV[Face[i][2]],dsize1);
		PT->setParametricDt(bc[1],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pU[Face[i][0]],pU[Face[i][1]],
				pU[Face[i][2]],dsize1);
    
		E[i] = PT->InnerProduct(bc[0],bc[0]);
		G[i] = PT->InnerProduct(bc[1],bc[1]);
	}
	PointTool pt;
	for(i=0;i<numberV;i++)
	{
         
		sigma[i]=0.0;
		now = FHead[i];
		varphi=0.0;
		localsum=0.0;	
		double varphiU = 0;
		double varphiV = 0;
		VList *nowv = VHead[i];
		Point3d *PointNormal = &(VHead[i]->normal);
		while (next(nowv) != VTail[i])
		{
			nowv = next(nowv);
			double x = E[nowv->FaceID]/G[nowv->FaceID];
			x = (x + (1.0/x));		
			double m = 2;
			double l = -(1.0/pow(2,m))*pow(x,m) + 1;
			double weight = exp(l);				
			

			 
			varphi  += (mirrorFace[nowv->FaceID].area3d *(0.5*(E[nowv->FaceID]+G[nowv->FaceID])));
			varphiV += (mirrorFace[nowv->FaceID].area3d *G[nowv->FaceID]);
			varphiU += (mirrorFace[nowv->FaceID].area3d *E[nowv->FaceID]);
			localsum += mirrorFace[nowv->FaceID].area3d; 
			nowv = next(nowv);
		}

		sigma[i] = sqrt((varphi/localsum));		
		sigma[i] = pow(sigma[i],gammaP);
		sigmaU[i] = sqrt((varphiU/localsum));
		sigmaU[i] = pow(sigmaU[i],gammaP);
		sigmaV[i] = sqrt((varphiV/localsum));
		sigmaV[i] = pow(sigmaV[i],gammaP);
	}  
}


void	MyParameterization::setSigmaU()
{
	int i;
	
	double dsize1=0.0;
	double dddhval=0.0;
	double localsum=0.0;
	for(i=0;i<numberF;i++)
	{
		dsize1 = PT->getParametricA(pV[Face[i][0]],
					pV[Face[i][1]],
					pV[Face[i][2]],
					pU[Face[i][0]],
					pU[Face[i][1]],
					pU[Face[i][2]]);
		
		PT->setParametricDs(bc[0],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pV[Face[i][0]],pV[Face[i][1]],
				pV[Face[i][2]],dsize1);
		E[i] = PT->InnerProduct(bc[0],bc[0]);
	}
	
	for(i=0;i<numberV;i++)
	{       
		localsum=0.0;	
		double varphiU = 0;		
		VList *nowv = VHead[i];
		while (next(nowv) != VTail[i])
		{
			nowv = next(nowv);
			varphiU += (mirrorFace[nowv->FaceID].area3d *E[nowv->FaceID]);
			localsum += mirrorFace[nowv->FaceID].area3d; 
			nowv = next(nowv);
		}	
		sigmaU[i] = sqrt((varphiU/localsum));
		sigmaU[i] = pow(sigmaU[i],gammaU);
	} 
}

void	MyParameterization::setSigmaV()
{
	int i;
	
	double dsize1=0.0;
	double dddhval=0.0;
	double localsum=0.0;
	for(i=0;i<numberF;i++)
	{
		dsize1 = PT->getParametricA(pV[Face[i][0]],
					pV[Face[i][1]],
					pV[Face[i][2]],
					pU[Face[i][0]],
					pU[Face[i][1]],
					pU[Face[i][2]]);		
		PT->setParametricDt(bc[1],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pU[Face[i][0]],pU[Face[i][1]],
				pU[Face[i][2]],dsize1);

		G[i] = PT->InnerProduct(bc[1],bc[1]);
	}
	
	for(i=0;i<numberV;i++)
	{       
		localsum=0.0;	
		double varphiV = 0;		
		VList *nowv = VHead[i];
		while (next(nowv) != VTail[i])
		{
			nowv = next(nowv);
			varphiV += (mirrorFace[nowv->FaceID].area3d *G[nowv->FaceID]);
			localsum += mirrorFace[nowv->FaceID].area3d; 
			nowv = next(nowv);
		}	
		sigmaV[i] = sqrt((varphiV/localsum));
		sigmaV[i] = pow(sigmaV[i],gammaV);
	} 
}


double MyParameterization::getCurrentE_EX(){
  int i,j;
  IDList *now=NULL;
  double varphi,ddv,dsize1,sumarea;
  double dddhval=0.0;
  double dsum=0.0;
  double dsumU=0.0;
  double dsumV=0.0;
  
  
  
  double localsum=0.0;
  sumarea = 0.0;
  double dE,dG;
  dE = 0.0;
  dG = 0.0;

  double max_strect = 0.0;
  if(boundarysigma==0)
  {
	  for(i=0;i<numberF;i++)
	  {
		if(boundary[Face[i][0]]!=1&&boundary[Face[i][1]]!=1&&boundary[Face[i][2]]!=1)
		{
			dsize1 = PT->getParametricA(pV[Face[i][0]],
					pV[Face[i][1]],
					pV[Face[i][2]],
					pU[Face[i][0]],
					pU[Face[i][1]],
					pU[Face[i][2]]);
		
			PT->setParametricDs(bc[0],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pV[Face[i][0]],pV[Face[i][1]],
				pV[Face[i][2]],dsize1);
			PT->setParametricDt(bc[1],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pU[Face[i][0]],pU[Face[i][1]],
				pU[Face[i][2]],dsize1);
    
			dE = PT->InnerProduct(bc[0],bc[0]);
			dG = PT->InnerProduct(bc[1],bc[1]);
			dsum +=  mirrorFace[i].area3d*0.5*(dE+dG);
			dsumV +=  mirrorFace[i].area3d*dE;
			dsumU +=  mirrorFace[i].area3d*dG;			
		}
	  }
	  sumarea += mirrorFace[i].area3d;	
	}
	else
	{
		for(i=0;i<numberF;i++)
		{
			dsize1 = PT->getParametricA(pV[Face[i][0]],
					pV[Face[i][1]],
					pV[Face[i][2]],
					pU[Face[i][0]],
					pU[Face[i][1]],
					pU[Face[i][2]]);
		
			PT->setParametricDs(bc[0],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pV[Face[i][0]],pV[Face[i][1]],
				pV[Face[i][2]],dsize1);
			PT->setParametricDt(bc[1],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pU[Face[i][0]],pU[Face[i][1]],
				pU[Face[i][2]],dsize1);
    
			dE = PT->InnerProduct(bc[0],bc[0]);
			dG = PT->InnerProduct(bc[1],bc[1]);
			dsum +=  mirrorFace[i].area3d*0.5*(dE+dG);
			dsumV +=  mirrorFace[i].area3d*dE;
			dsumU +=  mirrorFace[i].area3d*dG;
			sumarea += mirrorFace[i].area3d;
		}  
	}
	//dsum =constsumarea3D*sqrt(dsum/sumarea3D);  
	double constant_mirrorFaceArea = 0;
	if(boundarytype==0||boundarytype==2)
		constant_mirrorFaceArea = sqrt(1.0/sumarea);
	else if(boundarytype==1)	
		constant_mirrorFaceArea = sqrt(((0.5*0.5*PI)/sumarea));

	dsum = constant_mirrorFaceArea*sqrt(dsumU/sumarea) + constant_mirrorFaceArea*sqrt(dsumV/sumarea) ;
	return dsum;

}


double	MyParameterization::getCurrentE_U()
{
	int i,j;
	IDList *now=NULL;
	double dsize1,sumarea;
	
	
	double dsumU=0.0;
	double dE;
	dE = 0.0;
	
  
  
	double localsum=0.0;
	sumarea = 0.0;
	
	

  
	if(boundarysigma==0)
	{
		for(i=0;i<numberF;i++)
		{
		if(boundary[Face[i][0]]!=1&&boundary[Face[i][1]]!=1&&boundary[Face[i][2]]!=1)
		{
			dsize1 = PT->getParametricA(pV[Face[i][0]],
					pV[Face[i][1]],
					pV[Face[i][2]],
					pU[Face[i][0]],
					pU[Face[i][1]],
					pU[Face[i][2]]);
		
			PT->setParametricDs(bc[0],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pV[Face[i][0]],pV[Face[i][1]],
				pV[Face[i][2]],dsize1);
    
			dE = PT->InnerProduct(bc[0],bc[0]);			
			dsumU +=  mirrorFace[i].area3d*dE;			
		}
		}
		sumarea += mirrorFace[i].area3d;	
	}
	else
	{
		for(i=0;i<numberF;i++)
		{
			dsize1 = PT->getParametricA(pV[Face[i][0]],
					pV[Face[i][1]],
					pV[Face[i][2]],
					pU[Face[i][0]],
					pU[Face[i][1]],
					pU[Face[i][2]]);
		
			PT->setParametricDs(bc[0],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pV[Face[i][0]],pV[Face[i][1]],
				pV[Face[i][2]],dsize1);
    
			dE = PT->InnerProduct(bc[0],bc[0]);			
			dsumU +=  mirrorFace[i].area3d*dE;	
			sumarea += mirrorFace[i].area3d;
		}  
	}
	//dsum =constsumarea3D*sqrt(dsum/sumarea3D);  
	double constant_mirrorFaceArea = 0;
	if(boundarytype==0||boundarytype==2)
		constant_mirrorFaceArea = sqrt(1.0/sumarea);
	else if(boundarytype==1)	
		constant_mirrorFaceArea = sqrt(((0.5*0.5*PI)/sumarea));

	return constant_mirrorFaceArea*sqrt(dsumU/sumarea);	 
}

double	MyParameterization::getCurrentE_V()
{
	int i,j;
	IDList *now=NULL;
	double dsize1,sumarea;
	
	
	double dsumV=0.0;
	double dG;	
	dG = 0.0;
  
  
	double localsum=0.0;
	sumarea = 0.0;
	
	

  
	if(boundarysigma==0)
	{
		for(i=0;i<numberF;i++)
		{
		if(boundary[Face[i][0]]!=1&&boundary[Face[i][1]]!=1&&boundary[Face[i][2]]!=1)
		{
			dsize1 = PT->getParametricA(pV[Face[i][0]],
					pV[Face[i][1]],
					pV[Face[i][2]],
					pU[Face[i][0]],
					pU[Face[i][1]],
					pU[Face[i][2]]);
		
			PT->setParametricDt(bc[1],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pU[Face[i][0]],pU[Face[i][1]],
				pU[Face[i][2]],dsize1);
    
			dG = PT->InnerProduct(bc[1],bc[1]);			
			dsumV +=  mirrorFace[i].area3d*dG;
		}
		}
		sumarea += mirrorFace[i].area3d;	
	}
	else
	{
		for(i=0;i<numberF;i++)
		{
			dsize1 = PT->getParametricA(pV[Face[i][0]],
					pV[Face[i][1]],
					pV[Face[i][2]],
					pU[Face[i][0]],
					pU[Face[i][1]],
					pU[Face[i][2]]);
		
			PT->setParametricDt(bc[1],&mirrorFace[i].pt[0],
				&mirrorFace[i].pt[1],&mirrorFace[i].pt[2],
				pU[Face[i][0]],pU[Face[i][1]],
				pU[Face[i][2]],dsize1);
    
			dG = PT->InnerProduct(bc[1],bc[1]);			
			dsumV +=  mirrorFace[i].area3d*dG;	
			sumarea += mirrorFace[i].area3d;
		}  
	}
	//dsum =constsumarea3D*sqrt(dsum/sumarea3D);  
	double constant_mirrorFaceArea = 0;
	if(boundarytype==0||boundarytype==2)
		constant_mirrorFaceArea = sqrt(1.0/sumarea);
	else if(boundarytype==1)	
		constant_mirrorFaceArea = sqrt(((0.5*0.5*PI)/sumarea));

	return constant_mirrorFaceArea*sqrt(dsumV/sumarea);
}


void MyParameterization::ParametrizationOptimalSaveU0(double *op_u0, double *op_v0,int itenum,double error,FILE* logFile)
{

    double *UaXY = new double[2*(numberV)+1];
    double *vecb = new double[2*(numberV)+1];
  
	int i;
	IDList *now;
	IDList *now2;
	PolarList *nowp;
	level=0;
  
	int nonzero=(numberV);
	for(i=0;i<numberV;i++)
	{
		vecb[i+1]=0.0;
		if(boundary[i]!=1)
		{
			nonzero += neighborI[i];			
		}
	}

	int iter=0;
	double linerr=0.0;
	double weight=0.0;
  
	PCBCGSolver *mybcg = new PCBCGSolver(2*nonzero);
	double *sigsum = new double[numberV];
	if(weighttype==0)
	{
		setFloaterC();
	}
	else if(weighttype==1)
	{
		setLaplaceC();
	}
	else if(weighttype==2)
	{
		setEckHC();
	}
	else if(weighttype==3)
	{
		setDesbrunC();
	}
	else if(weighttype==4)
	{
		setMVCC();
	}
	else
	{
		setFloaterC();
	}
    
  
	SortIndexP();
  
  

	for(i=0;i<numberV;i++)
	{    
		if(boundary[i]!=1)
		{
			mybcg->sa[i+1] = 1.0;
			mybcg->sa[i+1+numberV] = 1.0;
			vecb[i+1] = 0.0;
			vecb[i+1+numberV] = 0.0;
		}
		else
		{
			mybcg->sa[i+1] = 1.0;
			vecb[i+1] = pU[i];
			mybcg->sa[i+1+numberV] = 1.0;
			vecb[i+1+numberV] = pV[i];
		}
	}

	mybcg->ija[1] = 2*(numberV)+2;
	int dlk=2*(numberV)+1;
  
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			nowp = PHead[i];      
			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);
				++dlk;
				mybcg->sa[dlk] = -nowp->lambda;
				mybcg->ija[dlk]=nowp->ID+1;
			}
		}
		mybcg->ija[i+1+1]=dlk+1;
	}
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			nowp = PHead[i];
			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);
				++dlk;
				mybcg->sa[dlk] = -nowp->lambda;
				mybcg->ija[dlk]=nowp->ID+numberV+1;
			}
		}
		mybcg->ija[i+numberV+1+1]=dlk+1;
	}
  
	for(i=0;i<numberV;i++)
	{
		UaXY[i+1] = pU[i];
		UaXY[i+numberV+1] = pV[i];
	}
	mybcg->linbcg(((unsigned long)(2*(numberV))),vecb,UaXY,1,error,itenum,&iter,&linerr);
	if (iter < itenum && linerr <= error )
	{
		for(i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				pU[i] = UaXY[i+1];
				pV[i] = UaXY[i+numberV+1];

				if (op_u0)
				{
					op_u0[i] = UaXY[i+1];
				}
				if (op_v0)
				{
					op_v0[i] = UaXY[i+numberV+1];
				}
			}
			/*
			else
			{
				if (op_u0)
				{
					op_u0[i] = pU[i];
				}

				if (op_v0)
				{
					op_v0[i] = pV[i];
				}
			}
			*/
		}
	}
	else
	{
		for(i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				pU[i] = 0.5;
				pV[i] = 0.5;
				if (op_u0)
				{
					op_u0[i] = 0.5;
				}
				if (op_v0)
				{
					op_v0[i] = 0.5;
				}
			}
		}
		resultStretch = DBL_MAX;
		printf("Fail to solve\n");
		delete mybcg;
		delete [] vecb;
		delete [] UaXY;
		delete [] sigsum;
		return;
	}
	// Re-solving linear system
	double initialstrech=0.0;
	double currentstrech=0.0;
	initialstrech = getCurrentE();
	
	if (initialstrech < 1)
	{
		resultStretch = DBL_MAX;
		printf("Fail to solve\n");
		delete mybcg;
		delete [] vecb;
		delete [] UaXY;
		delete [] sigsum;
		return;
	}

	
	
	int kk = 0;  
	//printf("U%d  STRETCH: %f\n",kk,initialstrech); //u0
	if (logFile)
	{
		fprintf(logFile,"U%d  STRETCH: %f\n",kk,initialstrech); //u0
	}

	Point2d *prevUV = new Point2d[numberV];	


	for(kk=0;kk<itenum;kk++)
	{
		setSigmaZero();
		for(i=0;i<numberV;i++)
		{
      
			if(boundary[i]!=1)
			{
				sigsum[i]=0.0;
				nowp = PHead[i];
				while(next(nowp)!=PTail[i])
				{
					nowp = next(nowp);
					nowp->old_lambda = nowp->lambda;
					nowp->lambda /= sigma[nowp->ID];
					sigsum[i] += nowp->lambda;
				}
				mybcg->sa[i+1] = 1.0;
				mybcg->sa[i+1+numberV] = 1.0;
			}
			else
			{
				mybcg->sa[i+1] = 1.0;
				mybcg->sa[i+1+numberV] = 1.0;
			}
		}
		dlk=2*(numberV)+1;
    
		for(i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				nowp = PHead[i];
	
				while(next(nowp)!=PTail[i])
				{
					nowp = next(nowp);
					++dlk;
					mybcg->sa[dlk] = -nowp->lambda/sigsum[i];
				}
			}
		}

		for(i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				nowp = PHead[i];
				while(next(nowp)!=PTail[i])
				{
					nowp = next(nowp);
					++dlk;
					mybcg->sa[dlk] = -nowp->lambda/sigsum[i];
				}
			}
		}
     
		for(i=0;i<numberV;i++)
		{
			UaXY[i+1] = pU[i];
			UaXY[i+numberV+1] = pV[i];
		}

		mybcg->linbcg(((unsigned long)(2*(numberV))),vecb,UaXY,1,error,itenum,&iter,&linerr); 
		if (iter < itenum)
		{
			for(i=0;i<numberV;i++)
			{
				if(boundary[i]!=1)
				{
					prevUV[i].x = pU[i];
					prevUV[i].y = pV[i];
	
					pU[i] = UaXY[i+1];
					pV[i] = UaXY[i+numberV+1];
				}
			}
			currentstrech = getCurrentE();	
		}
		else
		{
			currentstrech = -1;
		}
		
		//printf("currentstrech U%d= %lf\n",kk+1,currentstrech);
		//printf("U%d  STRETCH: %f\n",kk+1,currentstrech); 

		if (logFile)
		{
			fprintf(logFile,"U%d  STRETCH: %f\n",kk+1,currentstrech); 
		}

		if(initialstrech<currentstrech || currentstrech < 0)
		{
			for(i=0;i<numberV;i++)
			{
				if(boundary[i]!=1)
				{
					pU[i] = prevUV[i].x;
					pV[i] = prevUV[i].y;
					
				}
			}

			resultStretch = initialstrech;
			break;		
		}
		else
		{
			initialstrech = currentstrech;
		}
	}  
	level = kk;
 
 
	printf("STL2 error (U%d) = %lf\n",level,resultStretch);
	delete mybcg;
	
	delete [] prevUV;
	delete [] vecb;
	delete [] UaXY;
	delete [] sigsum;

}


/*
double    MyParameterization::PARAM_VIENNA_CL(	PolarVertex *pIPV,
												int num_PV,
												FILE* logFile)
{
	

	double bestStretch = -1;
	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;	
	if (pU)
		delete [] pU;
	if (pV)
		delete [] pV;

	pU = new double[numberV];
	pV = new double[numberV];
	memset(pU,0x00,sizeof(double)*numberV);
	memset(pV,0x00,sizeof(double)*numberV);
	
	double *initU0 = new double[numberV];
	double *initV0 = new double[numberV];
	memset(initU0,0x00,sizeof(double)*numberV);
	memset(initV0,0x00,sizeof(double)*numberV);

	clock_t start,finish;
	start = clock();

    boundarytype=0;
	setPolarMap();
	IDList *BpointH = new IDList();
	IDList *BpointT = new IDList();
    
	BpointH->next = BpointT;
	BpointT->back = BpointH;
	
	double tlength=0.0;
	int numBorderPoint = 0;
	CalBorderPath(BpointH,BpointT,&tlength,&numBorderPoint);

	IDList *now = BpointH;
	
	int state=0;
	double clen=0.0;		
	double loop = 0;
	double sum_length = 0;
	double edge_length = 0;		
	IDList *startPoint = BpointH->next;
	int best_startID = startPoint->ID;	
	int worst_startID = startPoint->ID;
	double worstStretch = -1;
	int count = 0;

	FILE *fpCorner = fopen("cornerID_setting.txt","rt");
	int spec_cornerID = -1;

	if (fpCorner)
	{
		printf("FOUND CORNER SETTING FILE!\n");
		//read setting from file		
		fscanf(fpCorner,"%d\r\n",&spec_cornerID);		
		fclose(fpCorner);
		if (spec_cornerID >= 0)
			printf("Calculate CORNER ID%d\n",spec_cornerID);
		else
		{
			spec_cornerID = -1;
			printf("Calculate ALL corner\n");
		}
		 
	}

	//record result as array.
	double *stretch_each_corner = new double[numBorderPoint];
	while (loop < tlength*0.25)
	{
		BpointT->ID = BpointH->next->ID; // for cal length;
		state = 0;
		sum_length = 0;
		double this_side_length[4] ={0};
		int this_side_num_edge[4] = {0};		
		loop += PT->Distance(point[startPoint->ID],point[next(startPoint)->ID]);
		IDList *now = startPoint;
		if (spec_cornerID < 0 || spec_cornerID == count)
		{
			do
			{				
				edge_length = PT->Distance(point[now->ID],point[next(now)->ID]);			
				if ((sum_length + edge_length>= (tlength*0.25)) && state < 3)
				{
					if ((sum_length + edge_length)/(tlength*0.25) >= 1.1&& this_side_num_edge[state] > 0)					
					{
						this_side_length[state] = sum_length;
						this_side_num_edge[state+1]++;
						sum_length = edge_length;
					}
					else
					{
						this_side_length[state] = sum_length+edge_length;
						this_side_num_edge[state]++;
						sum_length = 0;
					}				
					state++;
				}
				else
				{
					this_side_num_edge[state]++;
					sum_length += edge_length;
				}
				now = next(now);
				if (now == BpointT)
					now = BpointH->next;
			}
			while (now != startPoint);			
			this_side_length[state] = sum_length;

			//memcpy(pU,initU0,sizeof(double)*numberV);
			//memcpy(pV,initV0,sizeof(double)*numberV);

			state = 0;
			sum_length  = 0;
			clen = 0.0;
			now = BpointH;
			BpointT->ID = BpointH->next->ID; // for cal length;
		
			printf("TEST %d (%f/%f)  POINT CORNER ID: ",count,loop,tlength*0.25);		
			printf("%d ",BpointH->next->ID);
			
			if (logFile)
			{
				fprintf(logFile,"TEST %d (%f/%f)  POINT CORNER ID: ",count,loop,tlength*0.25);
				fprintf(logFile,"%d ",BpointH->next->ID);
			}
			
			memset(pU,0x00,sizeof(double)*numberV);
			memset(pV,0x00,sizeof(double)*numberV);	
			while(next(now)!=BpointT)
			{
				now = next(now);
				switch(state)
				{
					case 0:
						if (clen < 1.0)
						{
							pU[now->ID] = clen; 
							pV[now->ID] = 0.0;
						}
						else
						{
							pU[now->ID] = 1.0;
							pV[now->ID] = 0.0;
							sum_length = 0.0;
							state=1;
							printf("%d ",now->ID);
							
							if (logFile)
								fprintf(logFile,"%d ",now->ID);
							
						}
						break;
					case 1:
						if (clen < 1.0)
						{
							pU[now->ID] = 1.0; 
							pV[now->ID] = clen;
						}
						else
						{
							pU[now->ID] = 1.0;
							pV[now->ID] = 1.0;
							sum_length = 0.0;
							state=2;
							printf("%d ",now->ID);
							
							if (logFile)
								fprintf(logFile,"%d ",now->ID);
							

						}
						break;
					case 2:
						if (clen < 1.0)
						{
							pU[now->ID] = 1.0-clen; 
							pV[now->ID] = 1.0;
						}
						else
						{
							pU[now->ID] = 0.0;
							pV[now->ID] = 1.0;
							sum_length = 0.0;
							state=3;
							printf("%d\n",now->ID);
							
							if (logFile)
								fprintf(logFile,"%d\n",now->ID);
							

						}
						break;
					case 3:
						if (clen < 1.0)
						{
							pU[now->ID] = 0.0; 
							pV[now->ID] = 1.0-clen;
						}
						else
						{
							//should never enter this one;
							state=4;
						}
					default:
						break;

				}		
				sum_length += PT->Distance(point[now->ID],point[next(now)->ID]);
				clen = sum_length/this_side_length[state];	
			}
		
			
			ResetInnerLambda();		
			gammaP = 0.75;
			
			//ParametrizationOptimalSaveU0(initU0,initV0,iteNum,PCBCGerror,NULL);

			this->ParametrizationOptimalGPU(PCBCGerror,NULL);
			printf("STRETCH=%f\n",resultStretch);
			if (logFile)
			{
				fprintf(logFile,"STRETCH=%f\n",resultStretch);
				
			}
			gammaP = 1.0;
			if (resultStretch < bestStretch || bestStretch < 0)
			{
				for (int i=0;i<numberV;i++)
				{
					if(boundary[i]!=1)
					{
						PolarList *nowp = PHead[i];
						while(next(nowp)!=PTail[i])
						{
							nowp = next(nowp);
							nowp->best_lambda = nowp->old_lambda;
						}
					}
					pIPV[i].u = pU[i];
					pIPV[i].v = pV[i];
				}
				bestStretch = resultStretch;
				best_startID = count;
			}

			if (resultStretch > worstStretch)
			{
				worstStretch = resultStretch;
				worst_startID = count;
			}
		}

		if (loop >= tlength*0.25)
			break;

		count++;
		//prepare for next loop
		startPoint = next(startPoint);
		//reconstruct BpointH , BpointT according to startPoint
		int startPointID = startPoint->ID;
		while (BpointH->next->ID != startPoint->ID)
		{
			IDList *move_node = BpointH->next;
			BpointH->next = move_node->next;
			move_node->next->back = BpointH;
			
			move_node->back = BpointT->back;
			BpointT->back->next = move_node;
			BpointT->back = move_node;
			move_node->next = BpointT;
		}
	}
	
	delete [] initU0;
	delete [] initV0;


	finish = clock();
	printf( "Calculate Time for best corner : %f sec\n", (double)(finish-start)/CLOCKS_PER_SEC);
	if (logFile)
		fprintf(logFile,"Calculate Time for best corner : %f sec\n", (double)(finish-start)/CLOCKS_PER_SEC);
	

	for (int i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			PolarList *nowp = PHead[i];
			while(next(nowp)!=PTail[i])
			{
				nowp = next(nowp);
				nowp->lambda = nowp->best_lambda;
			}
		}
		pU[i] = pIPV[i].u;
		pV[i] = pIPV[i].v;
	}
	
	printf("=== Find best stretch  of TEST%d  (best corner at %f) ===\n",best_startID,bestStretch);
	if (logFile)
		fprintf(logFile,"=== Find best stretch of TEST%d (best corner at %f)===\n",best_startID,bestStretch);

	ParametrizationSmoothOptimal_EX(0.5,0.1,iteNum,PCBCGerror,logFile);

	printf("=== Find worst stretch  of TEST%d  (worst corner at %f) ===\n",worst_startID,worstStretch);
	if (logFile)
		fprintf(logFile,"=== Find worst stretch  of TEST%d  (worst corner at %f) ===\n",worst_startID,worstStretch);
	
	for (int i=0;i<numberV;i++)
	{
		pIPV[i].u = pU[i];
		pIPV[i].v = pV[i];
	}
	IDtool->CleanNeighbor(BpointH,BpointT);
	return resultStretch;



	
}
void MyParameterization::SquareParametrizationGPU(int bottomleftIDvertex,IDList *borderH,IDList *borderT,double total_length_edges, cpuCompressedMatrixType *initAMatrix,
											      double *& resultU,double *& resultV,double &resultStretchErr,bool directsolver,bool reportlog )
{
	vector<int> bordervertexT;
	vector<int> bordervertex;
	vector<int>::iterator bv_iter;
	IDList *now=borderH;
	int count = 0;
	while(next(now)!=borderT)
	{
		now = next(now);
		if (now->ID == bottomleftIDvertex)
		{
			bordervertex.push_back(now->ID);
			while(next(now)!=borderT)
			{
				now = next(now);
				bordervertex.push_back(now->ID);
			}
		}
		else
		{
			bordervertexT.push_back(now->ID);
			count++;
		}
	}	
	if (bordervertexT.size()>0)
		bordervertex.insert(bordervertex.end(), bordervertexT.begin(), bordervertexT.end());
	
	bordervertex.push_back(bottomleftIDvertex); //for easy loop back at end

	
	int state = 0;
	double this_side_length[4] ={0};
	int this_side_num_edge[4] = {0};
	double sum_length = 0;
	bv_iter = bordervertex.begin();
	do
	{				
		double edge_length = PT->Distance(point[*bv_iter],point[*(bv_iter+1)]);			
		if ((sum_length + edge_length>= (total_length_edges*0.25)) && state < 3)
		{
			if ((sum_length + edge_length)/(total_length_edges*0.25) >= 1.1&& this_side_num_edge[state] > 0)					
			{
				this_side_length[state] = sum_length;
				this_side_num_edge[state+1]++;
				sum_length = edge_length;
			}
			else
			{
				this_side_length[state] = sum_length+edge_length;
				this_side_num_edge[state]++;
				sum_length = 0;
			}				
			state++;
		}
		else
		{
			this_side_num_edge[state]++;
			sum_length += edge_length;
		}
		bv_iter++;
	}
	while ((bv_iter) != bordervertex.end() - 1);
	this_side_length[state] = sum_length;  //for last stage
	
	//memcpy(pU,initU0,sizeof(double)*numberV);
	//memcpy(pV,initV0,sizeof(double)*numberV);

	state = 0;
	sum_length  = 0;
	double clen = 0.0;
	if (reportlog)
	{
		printf("TEST %d POINT CORNER ID: ",count);		
		printf("%d ",bordervertex[0]);
	}
			
	if (logFile)
	{
		fprintf(logFile,"TEST %d (%f/%f)  POINT CORNER ID: ",count,loop,tlength*0.25);
		fprintf(logFile,"%d ",BpointH->next->ID);
	}
	
	
	double *ioU = new double [numberV];
	double *ioV = new double [numberV];
	memset(ioU,0x00,sizeof(double)*numberV);
	memset(ioV,0x00,sizeof(double)*numberV);	
	bv_iter = bordervertex.begin();
	do
	{	
		switch(state)
		{
			case 0:
				if (clen < 1.0)
				{
					ioU[*bv_iter] = clen; 
					ioV[*bv_iter] = 0.0;
				}
				else
				{
					ioU[*bv_iter] = 1.0;
					ioV[*bv_iter] = 0.0;
					sum_length = 0.0;
					state=1;
					if (reportlog)	
						printf("%d ",*bv_iter);
							
					//if (logFile)
					//	fprintf(logFile,"%d ",now->ID);
						
				}
				break;
			case 1:
				if (clen < 1.0)
				{
					ioU[*bv_iter] = 1.0; 
					ioV[*bv_iter] = clen;
				}
				else
				{
					ioU[*bv_iter] = 1.0;
					ioV[*bv_iter] = 1.0;
					sum_length = 0.0;
					state=2;
					if (reportlog)	
						printf("%d ",*bv_iter);
							
					//if (logFile)
					//	fprintf(logFile,"%d ",now->ID);
							

				}
				break;
			case 2:
				if (clen < 1.0)
				{
					ioU[*bv_iter] = 1.0-clen; 
					ioV[*bv_iter] = 1.0;
				}
				else
				{
					ioU[*bv_iter] = 0.0;
					ioV[*bv_iter] = 1.0;
					sum_length = 0.0;
					state=3;
					if (reportlog)
						printf("%d\n",*bv_iter);
							
					//if (logFile)
					//	fprintf(logFile,"%d\n",now->ID);
							

				}
				break;
			case 3:
				if (clen < 1.0)
				{
					ioU[*bv_iter] = 0.0; 
					ioV[*bv_iter] = 1.0-clen;
				}
				else
				{	
					state=4;
				}
			default:
				break;
		}
		sum_length += PT->Distance(point[*bv_iter],point[*(bv_iter+1)]);
		clen = sum_length/this_side_length[state];	
		bv_iter++;
	}
	while (bv_iter != bordervertex.end()-1);


	double stretchError = ParametrizationOptimalGPU( ioU,ioV,PCBCGerror,initAMatrix,directsolver,NULL);
	if (reportlog)	
		printf("stretchError = %f\n",stretchError);
	resultU = ioU;
	resultV = ioV;
	resultStretchErr = stretchError;

}
*/
void MyParameterization::ParametrizationSmoothOptimal_EX(double startGamma, double finishGamma,int itenum,double error,FILE* logFile)
{
	
	double *UaXY = new double[2*(numberV)+1];
	double *vecb = new double[2*(numberV)+1];
  
	int i;
	IDList *now;
	IDList *now2;
	PolarList *nowp;
	level=0;
  
	int nonzero=(numberV);
	for(i=0;i<numberV;i++)
	{		
		if(boundary[i]!=1)
		{
			nonzero += neighborI[i];
			vecb[i+1]=0.0;
			vecb[i+1+numberV] = 0.0;
		}
		else
		{
			vecb[i+1] = pU[i];			
			vecb[i+1+numberV] = pV[i];
		}
	}
	int iter=0;
	double linerr=0.0;
	double weight=0.0;
  
	PCBCGSolver *mybcg = new PCBCGSolver(2*nonzero);
	double *sigsum = new double[numberV]; 
	mybcg->ija[1] = 2*(numberV)+2;
	int dlk=2*(numberV)+1;  
  
	double initialstrech= getCurrentE();
	double currentstrech=0.0;
	Point2d *prevU = new Point2d [numberV];  
	initialstrech = getCurrentE();
 
	int kk = 0;
	printf("U%d  STRETCH: %f\n",kk,initialstrech); //u0
	if (logFile)
	{
		fprintf(logFile,"U%d  STRETCH: %f\n",kk,initialstrech); //u0
	}

	resultStretch = initialstrech;
	double gamma = startGamma;
	for(kk=0;kk<itenum;kk++)
	{		
		setSigma(gamma);
		for(i=0;i<numberV;i++)
		{      
			if(boundary[i]!=1)
			{
				sigsum[i]=0.0;
				nowp = PHead[i];
				while(next(nowp)!=PTail[i])
				{
					nowp = next(nowp);	  
					nowp->old_lambda = nowp->lambda;
					//nowp->lambda /= (0.5*(sigma[nowp->ID]+sigma[i]));
					nowp->lambda /= sigma[nowp->ID];
					sigsum[i] += nowp->lambda;
				}
				

				mybcg->sa[i+1] = 1.0;
				mybcg->sa[i+1+numberV] = 1.0;
			}	
			else
			{
				mybcg->sa[i+1] = 1.0;
				mybcg->sa[i+1+numberV] = 1.0;
			}
		}    
		dlk=2*(numberV)+1;
    
		for(i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				nowp = PHead[i];	
				while(next(nowp)!=PTail[i])
				{
					nowp = next(nowp);
					++dlk;
					mybcg->sa[dlk] = -nowp->lambda/sigsum[i];
					mybcg->ija[dlk]= nowp->ID+1;
				}
			}
			mybcg->ija[i+1+1]=dlk+1;
		}
		for(i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				nowp = PHead[i];
				while(next(nowp)!=PTail[i])
				{
					nowp = next(nowp);
					++dlk;
					mybcg->sa[dlk] = -nowp->lambda/sigsum[i];
					mybcg->ija[dlk]=nowp->ID+numberV+1 ;
				}
			}
			mybcg->ija[i+numberV+1+1]=dlk+1;
		}    
		for(i=0;i<numberV;i++)
		{
			if(boundary[i]==1)
			{
				UaXY[i+1] = pU[i];
				UaXY[i+numberV+1] = pV[i];
			}
			else
			{
				UaXY[i+1] = 0.5;
				UaXY[i+numberV+1] = 0.5;
			}
		}
		mybcg->linbcg(((unsigned long)(2*(numberV))),vecb,UaXY,1,error,itenum,&iter,&linerr);
		for(i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				prevU[i].x = pU[i];
				prevU[i].y = pV[i];
	
				pU[i] = UaXY[i+1];
				pV[i] = UaXY[i+numberV+1];
			}
		}
		currentstrech = getCurrentE();
		
		printf("U%d  STRETCH: %f\n",kk+1,currentstrech); //u1
		if (logFile)
		{
			fprintf(logFile,"U%d  STRETCH: %f\n",kk+1,currentstrech); //u1
		}
    
		if(initialstrech<currentstrech)
		{		
			for(i=0;i<numberV;i++)
			{
				if(boundary[i]!=1)
				{
					pU[i] = prevU[i].x;
					pV[i] = prevU[i].y;
				}
			}

			if (gamma > finishGamma)
			{
				gamma /= 2.0;
				if (gamma < finishGamma)
					gamma = finishGamma;
				for(i=0;i<numberV;i++)
				{      
					if(boundary[i]!=1)
					{
						nowp = PHead[i];
						while(next(nowp)!=PTail[i])
						{
							nowp = next(nowp);	  
							nowp->lambda = nowp->old_lambda;
						}
					}
				}    
			}
			else
			{			
				resultStretch = initialstrech;
				break;
			}			
		}
		else
		{
			initialstrech = currentstrech;
		}
	}
  
	level = kk;
  
  
	printf("U%d STL2 error = %lf\n",level,resultStretch);
	if (logFile)
	{
		fprintf(logFile,"U%d STL2 error = %lf\n",level,resultStretch);
	}
	


	delete mybcg;
	delete [] prevU;
	delete [] vecb;
	delete [] UaXY;
	delete [] sigsum;

}
void MyParameterization::SetSigma(double *ipU,double *ipV,double *opSigma,double gamma )
{
	IDList *now=NULL;
	double varphi,ddv,dsize1,sumarea;
	double dddhval=0.0;
	double localsum=0.0;
	Point3d _bc[2];
	double *_E = new double[numberF];
	double *_G = new double[numberF];
	for(int i=0;i<numberF;i++)
	{    
		dsize1 = PT->getParametricA(ipV[Face[i][0]],
									ipV[Face[i][1]],
									ipV[Face[i][2]],
									ipU[Face[i][0]],
									ipU[Face[i][1]],
									ipU[Face[i][2]]);
		PT->setParametricDs(&_bc[0],
							point[Face[i][0]],point[Face[i][1]],point[Face[i][2]],
							ipV[Face[i][0]],ipV[Face[i][1]],ipV[Face[i][2]],dsize1);
		PT->setParametricDt(&_bc[1],
							point[Face[i][0]],point[Face[i][1]],point[Face[i][2]],
							ipU[Face[i][0]],ipU[Face[i][1]],ipU[Face[i][2]],dsize1);
    
		_E[i] = PT->InnerProduct(&_bc[0],&_bc[0]);
		_G[i] = PT->InnerProduct(&_bc[1],&_bc[1]);
	}

	for(int i=0;i<numberV;i++)
	{
         
		opSigma[i]=0.0;
		now = FHead[i];
		varphi=0.0;
		localsum=0.0;   
    
		while(next(now)!=FTail[i])
		{
			now = next(now);
			varphi += (areaMap3D[now->ID]*(0.5*(_E[now->ID]+_G[now->ID])));
			localsum += (areaMap3D[now->ID]);
		}
    
		opSigma[i] = sqrt((varphi/localsum));    
		opSigma[i] = pow(opSigma[i],gamma);  
	}
	delete [] _E;
	delete [] _G;
  
}

void MyParameterization::setSigma(double gamma)
{
  int i,j;



#pragma omp parallel for
  for(i=0;i<numberF;i++){
    Point3d _bc[2];
    double dsize1 = PT->getParametricA(pV[Face[i][0]],
				pV[Face[i][1]],
				pV[Face[i][2]],
				pU[Face[i][0]],
				pU[Face[i][1]],
				pU[Face[i][2]]);
    PT->setParametricDs(&_bc[0],point[Face[i][0]],
			point[Face[i][1]],point[Face[i][2]],
			pV[Face[i][0]],pV[Face[i][1]],
			pV[Face[i][2]],dsize1);
    PT->setParametricDt(&_bc[1],point[Face[i][0]],
			point[Face[i][1]],point[Face[i][2]],
			pU[Face[i][0]],pU[Face[i][1]],
			pU[Face[i][2]],dsize1);
    
    E[i] = PT->InnerProduct(&_bc[0],&_bc[0]);
    G[i] = PT->InnerProduct(&_bc[1],&_bc[1]);
  }

#pragma omp parallel for
  for(i=0;i<numberV;i++){
         
    sigma[i]=0.0;
    IDList * now = FHead[i];
    double varphi=0.0;
    double localsum=0.0;   
    
    while(next(now)!=FTail[i]){
      now = next(now);
      varphi += (areaMap3D[now->ID]*(0.5*(E[now->ID]+G[now->ID])));
      localsum += (areaMap3D[now->ID]);
      
      
    }
    
    sigma[i] = sqrt((varphi/localsum));
    
    sigma[i] = pow(sigma[i],gamma);
    
  
  }
  
  
  
} 

double MyParameterization::GetStretchError(double *ipU,double *ipV, bool include_boundary, double *opFaceStretch)
{
	int i;
	double dsum = 0;
	if(!include_boundary)
	{
#pragma omp parallel for
		for(i=0;i<numberF;i++)
		{
			if(boundary[Face[i][0]]!=1&&boundary[Face[i][1]]!=1&&boundary[Face[i][2]]!=1)
			{
				Point3d _bc[2];
				double pV1 = ipV[Face[i][0]];
				double pV2 = ipV[Face[i][1]];
				double pV3 = ipV[Face[i][2]];
				double pU1 = ipU[Face[i][0]];
				double pU2 = ipU[Face[i][1]];
				double pU3 = ipU[Face[i][2]];    
				double dsize1 = PT->getParametricA(pV1,pV2,pV3,pU1,pU2,pU3);      
				PT->setParametricDs(&_bc[0],point[Face[i][0]],
						point[Face[i][1]],point[Face[i][2]],
						pV1,pV2,pV3,dsize1);
				PT->setParametricDt(&_bc[1],point[Face[i][0]],
						point[Face[i][1]],point[Face[i][2]],
						pU1,pU2,pU3,dsize1);
				double dE = PT->InnerProduct(&_bc[0],&_bc[0]);    
				double dG = PT->InnerProduct(&_bc[1],&_bc[1]);
				opFaceStretch[i] = sqrt((dE+dG)*0.5);
				#pragma omp critical
				{
					dsum +=  areaMap3D[i]*0.5*(dE+dG);
				}
			}
		}
	}
	else
	{
#pragma omp parallel for
		for(i=0;i<numberF;i++)
		{      
			Point3d _bc[2];
			double pV1 = ipV[Face[i][0]];
			double pV2 = ipV[Face[i][1]];
			double pV3 = ipV[Face[i][2]];
			double pU1 = ipU[Face[i][0]];
			double pU2 = ipU[Face[i][1]];
			double pU3 = ipU[Face[i][2]];    
			double dsize1 = PT->getParametricA(pV1,pV2,pV3,pU1,pU2,pU3);      
			PT->setParametricDs(&_bc[0],point[Face[i][0]],
					point[Face[i][1]],point[Face[i][2]],
					pV1,pV2,pV3,dsize1);
			PT->setParametricDt(&_bc[1],point[Face[i][0]],
					point[Face[i][1]],point[Face[i][2]],
					pU1,pU2,pU3,dsize1);
			double dE = PT->InnerProduct(&_bc[0],&_bc[0]);    
			double dG = PT->InnerProduct(&_bc[1],&_bc[1]);
			opFaceStretch[i] = sqrt((dE+dG)*0.5);
			#pragma omp critical
			{
				dsum +=  areaMap3D[i]*0.5*(dE+dG);
			}
		}  
	}
	dsum =constsumarea3D*sqrt(dsum/sumarea3D);
	if (dsum == dsum)
		return dsum;
	else
		return -1.0;
}
double MyParameterization::GetStretchError(double *ipU,double *ipV)
{
	int i;
	double dsum = 0;
	if(boundarysigma==0)
	{
#pragma omp parallel for
		for(i=0;i<numberF;i++)
		{
			if(boundary[Face[i][0]]!=1&&boundary[Face[i][1]]!=1&&boundary[Face[i][2]]!=1)
			{
				Point3d _bc[2];
				double pV1 = ipV[Face[i][0]];
				double pV2 = ipV[Face[i][1]];
				double pV3 = ipV[Face[i][2]];
				double pU1 = ipU[Face[i][0]];
				double pU2 = ipU[Face[i][1]];
				double pU3 = ipU[Face[i][2]];    
				double dsize1 = PT->getParametricA(pV1,pV2,pV3,pU1,pU2,pU3);      
				PT->setParametricDs(&_bc[0],point[Face[i][0]],
						point[Face[i][1]],point[Face[i][2]],
						pV1,pV2,pV3,dsize1);
				PT->setParametricDt(&_bc[1],point[Face[i][0]],
						point[Face[i][1]],point[Face[i][2]],
						pU1,pU2,pU3,dsize1);
				double dE = PT->InnerProduct(&_bc[0],&_bc[0]);    
				double dG = PT->InnerProduct(&_bc[1],&_bc[1]);
				#pragma omp critical
				{
					dsum +=  areaMap3D[i]*0.5*(dE+dG);
				}
			}
		}
	}
	else
	{
#pragma omp parallel for
		for(i=0;i<numberF;i++)
		{      
			Point3d _bc[2];
			double pV1 = ipV[Face[i][0]];
			double pV2 = ipV[Face[i][1]];
			double pV3 = ipV[Face[i][2]];
			double pU1 = ipU[Face[i][0]];
			double pU2 = ipU[Face[i][1]];
			double pU3 = ipU[Face[i][2]];    
			double dsize1 = PT->getParametricA(pV1,pV2,pV3,pU1,pU2,pU3);      
			PT->setParametricDs(&_bc[0],point[Face[i][0]],
					point[Face[i][1]],point[Face[i][2]],
					pV1,pV2,pV3,dsize1);
			PT->setParametricDt(&_bc[1],point[Face[i][0]],
					point[Face[i][1]],point[Face[i][2]],
					pU1,pU2,pU3,dsize1);
			double dE = PT->InnerProduct(&_bc[0],&_bc[0]);
    
			double dG = PT->InnerProduct(&_bc[1],&_bc[1]);
			#pragma omp critical
			{
				dsum +=  areaMap3D[i]*0.5*(dE+dG);
			}
		}  
	}
	dsum =constsumarea3D*sqrt(dsum/sumarea3D);
	if (dsum == dsum)
		return dsum;
	else
		return -1.0;
	
}

void MyParameterization::SquareParametrizationCPU(int bottomleftIDvertex,IDList *borderH,IDList *borderT,double total_length_edges, int non_zero_element,double *init_sa,unsigned long *init_ija , double *& resultU,double *& resultV,double &resultStretchErr,bool reportlog)
{
	vector<int> bordervertexT;
	vector<int> bordervertex;
	vector<int>::iterator bv_iter;
	IDList *now=borderH;
	int count = 0;
	while(next(now)!=borderT)
	{
		now = next(now);
		if (now->ID == bottomleftIDvertex)
		{
			bordervertex.push_back(now->ID);
			while(next(now)!=borderT)
			{
				now = next(now);
				bordervertex.push_back(now->ID);
			}
		}
		else
		{
			bordervertexT.push_back(now->ID);
			count++;
		}
	}	
	if (bordervertexT.size()>0)
		bordervertex.insert(bordervertex.end(), bordervertexT.begin(), bordervertexT.end());
	
	bordervertex.push_back(bottomleftIDvertex); //for easy loop back at end

	
	int state = 0;
	double this_side_length[4] ={0};
	int this_side_num_edge[4] = {0};
	double sum_length = 0;
	bv_iter = bordervertex.begin();
	do
	{				
		double edge_length = PT->Distance(point[*bv_iter],point[*(bv_iter+1)]);			
		if ((sum_length + edge_length>= (total_length_edges*0.25)) && state < 3)
		{
			if ((sum_length + edge_length)/(total_length_edges*0.25) >= 1.1&& this_side_num_edge[state] > 0)					
			{
				this_side_length[state] = sum_length;
				this_side_num_edge[state+1]++;
				sum_length = edge_length;
			}
			else
			{
				this_side_length[state] = sum_length+edge_length;
				this_side_num_edge[state]++;
				sum_length = 0;
			}				
			state++;
		}
		else
		{
			this_side_num_edge[state]++;
			sum_length += edge_length;
		}
		bv_iter++;
	}
	while ((bv_iter) != bordervertex.end() - 1);
	this_side_length[state] = sum_length;  //for last stage
	
	//memcpy(pU,initU0,sizeof(double)*numberV);
	//memcpy(pV,initV0,sizeof(double)*numberV);

	state = 0;
	sum_length  = 0;
	double clen = 0.0;
	if (reportlog)
	{
		printf("TEST %d POINT CORNER ID: ",count);		
		printf("%d ",bordervertex[0]);
	}
	/*		
	if (logFile)
	{
		fprintf(logFile,"TEST %d (%f/%f)  POINT CORNER ID: ",count,loop,tlength*0.25);
		fprintf(logFile,"%d ",BpointH->next->ID);
	}
	
	*/
	double *ioU = new double [numberV];
	double *ioV = new double [numberV];
	memset(ioU,0x00,sizeof(double)*numberV);
	memset(ioV,0x00,sizeof(double)*numberV);	
	bv_iter = bordervertex.begin();
	do
	{	
		switch(state)
		{
			case 0:
				if (clen < 1.0)
				{
					ioU[*bv_iter] = clen; 
					ioV[*bv_iter] = 0.0;
				}
				else
				{
					ioU[*bv_iter] = 1.0;
					ioV[*bv_iter] = 0.0;
					sum_length = 0.0;
					state=1;
					if (reportlog)	
						printf("%d ",*bv_iter);
							
					//if (logFile)
					//	fprintf(logFile,"%d ",now->ID);
						
				}
				break;
			case 1:
				if (clen < 1.0)
				{
					ioU[*bv_iter] = 1.0; 
					ioV[*bv_iter] = clen;
				}
				else
				{
					ioU[*bv_iter] = 1.0;
					ioV[*bv_iter] = 1.0;
					sum_length = 0.0;
					state=2;
					if (reportlog)	
						printf("%d ",*bv_iter);
							
					//if (logFile)
					//	fprintf(logFile,"%d ",now->ID);
							

				}
				break;
			case 2:
				if (clen < 1.0)
				{
					ioU[*bv_iter] = 1.0-clen; 
					ioV[*bv_iter] = 1.0;
				}
				else
				{
					ioU[*bv_iter] = 0.0;
					ioV[*bv_iter] = 1.0;
					sum_length = 0.0;
					state=3;
					if (reportlog)
						printf("%d\n",*bv_iter);
							
					//if (logFile)
					//	fprintf(logFile,"%d\n",now->ID);
							

				}
				break;
			case 3:
				if (clen < 1.0)
				{
					ioU[*bv_iter] = 0.0; 
					ioV[*bv_iter] = 1.0-clen;
				}
				else
				{	
					state=4;
				}
			default:
				break;
		}
		sum_length += PT->Distance(point[*bv_iter],point[*(bv_iter+1)]);
		clen = sum_length/this_side_length[state];	
		bv_iter++;
	}
	while (bv_iter != bordervertex.end()-1);


	double stretchError = ParametrizationOptimalCPU( ioU,ioV,PCBCGerror,non_zero_element,init_sa,init_ija,NULL);
	
	if (reportlog)	
		printf("stretchError = %f\n",stretchError);
	resultU = ioU;
	resultV = ioV;
	resultStretchErr = stretchError;
}

double  MyParameterization::PARAM_PARALLEL_CPU(	PolarVertex *pIPV,
												int num_PV,
												FILE* logFile)
{
	
	if (pU)
	{
		delete [] pU;
		pU = NULL;
	}
	if (pV)
	{
		delete [] pV;
		pV = NULL;
	}


    boundarytype=0;
	setPolarMap();
	IDList *BpointH = new IDList();
	IDList *BpointT = new IDList();
    
	BpointH->next = BpointT;
	BpointT->back = BpointH;
	
	double tlength=0.0;
	int numBorderPoint = 0;
	CalBorderPath(BpointH,BpointT,&tlength,&numBorderPoint);

	IDList *now = BpointH;
	


	double loop = 0;

	double bestStretch = DBL_MAX;
	int best_startID = -1;	
	int worst_startID = -1;
	double worstStretch = -1;
	


	//ResetInnerLambda();
	setFloaterC();
	SortIndexP();
	//record result as array.
	vector<int> bottomLeftVertexIDList;
	//double *stretch_each_corner = new double[numBorderPoint];
	printf("=== NUMBER BORDER EDGES : %d ===\n",numBorderPoint);
	if (logFile)
		fprintf(logFile,"=== NUMBER BORDER EDGES : %d ===\n",numBorderPoint);
	while (loop < tlength*0.25 && next(now) != BpointT)
	{
		now = next(now);		
		loop += PT->Distance(point[now->ID],point[next(now)->ID]);
		bottomLeftVertexIDList.push_back(now->ID);
	}

	//create initial A matrix  ,it is same for all condition
	//to boost up speed
	double *init_sa = NULL;
	unsigned long *init_ija = NULL;
	int nonzero=(numberV);
	for(int i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			nonzero += neighborI[i];			
		}
	}	
	init_ija = new unsigned long[2*nonzero+2];
    init_sa = new double[2*nonzero+2];	
	init_ija[1] = 2*(numberV)+2;
	int dlk=2*(numberV)+1;
  
	for(int i=0;i<numberV;i++)
	{
		init_sa[i+1] = 1.0;		
		if(boundary[i]!=1)
		{			
			PolarList *nowp = PHead[i];      
			while(nextN(nowp)!=PTail[i])
			{
				nowp = nextN(nowp);
				++dlk;
				init_sa[dlk] = -nowp->lambda;
				init_ija[dlk]= nowp->ID+1;
			}
		}
		init_ija[i+1+1]=dlk+1;
	}
	for(int i=0;i<numberV;i++)
	{
		init_sa[i+1+numberV] = 1.0;
		if(boundary[i]!=1)
		{
			PolarList *nowp = PHead[i];
			while(nextN(nowp)!=PTail[i])
			{
				nowp = nextN(nowp);
				++dlk;
				init_sa[dlk] = -nowp->lambda;
				init_ija[dlk]=nowp->ID+numberV+1;
			}
		}
		init_ija[i+numberV+1+1]=dlk+1;
	}
	//prepare for store result array
	double *bestU = NULL;
	double *bestV = NULL;
	std::vector<double *>_resultU(bottomLeftVertexIDList.size(),NULL);
	std::vector<double *>_resultV(bottomLeftVertexIDList.size(),NULL);
	std::vector<double>_resultError(bottomLeftVertexIDList.size(),0.0);
	
	unsigned int n = bottomLeftVertexIDList.size(); 
#if 0
	concurrency::SchedulerPolicy oldpolicy = concurrency::CurrentScheduler::GetPolicy();	
	concurrency::SchedulerPolicy policy(oldpolicy);
	if (policy.GetPolicyValue(concurrency::MaxConcurrency) > 10)
		policy.SetConcurrencyLimits(1,10);	
	
	concurrency::CurrentScheduler::Create(policy);	
	concurrency::parallel_for(0u, n, [&_resultU,&_resultV,&_resultError,m_nonzero,&init_sa,&init_ija,bottomLeftVertexIDList,BpointH,BpointT,tlength,this](unsigned int i)  
	{
		SquareParametrizationCPU(bottomLeftVertexIDList[i], BpointH,BpointT, tlength,nonzero,init_sa,init_ija,_resultU[i],_resultV[i],_resultError[i],false);		
	});
	concurrency::CurrentScheduler::Detach();
#else
#pragma omp parallel for
	for (int i = 0 ; i < bottomLeftVertexIDList.size(); i++)
	{		
		SquareParametrizationCPU(bottomLeftVertexIDList[i], BpointH,BpointT, tlength,nonzero,init_sa,init_ija,_resultU[i],_resultV[i],_resultError[i],false);	
	}
#endif

	for (int i = 0 ; i < bottomLeftVertexIDList.size(); i++)
	{
		if (_resultError[i] < bestStretch )
		{
			if (bestU)
				delete [] bestU;
			if (bestV)
				delete [] bestV;
			bestU = _resultU[i];
			bestV = _resultV[i];
			
			bestStretch = _resultError[i];
			best_startID = i;
		}
		else
		{
			
			delete [] _resultU[i];
			delete [] _resultV[i];
		}

		if (_resultError[i] > worstStretch)
		{
			worstStretch = _resultError[i];
			worst_startID = i;
		}	
	}
	
	for (int i=0;i<numberV;i++)
	{				
		pIPV[i].u = bestU[i];
		pIPV[i].v = bestV[i];
	}

	
	double *resultCircleU = NULL;
	double *resultCircleV = NULL;
	double resultCircleErr = 0;
	CircleParametrizationCPU(BpointH,BpointT, tlength,nonzero,init_sa,init_ija,resultCircleU,resultCircleV,resultCircleErr,false);


	if (bestU)
		delete [] bestU;
	if (bestV)
		delete [] bestV;

	if (resultCircleU)
		delete [] resultCircleU;
	if (resultCircleV)
		delete [] resultCircleV;

	if (init_ija)
		delete [] init_ija;
	if (init_sa)
		delete [] init_sa;

	
	printf("=== Find best stretch  of TEST%d  (best corner at %f) ===\n",best_startID,bestStretch);
	if (logFile)
		fprintf(logFile,"=== Find best stretch of TEST%d (best corner ERR = %f)===\n",best_startID,bestStretch);
	
	printf("=== Find worst stretch  of TEST%d  (worst corner at %f) ===\n",worst_startID,worstStretch);
	if (logFile)
		fprintf(logFile,"=== Find worst stretch  of TEST%d  (worst corner ERR = %f) ===\n",worst_startID,worstStretch);

	printf("=== Circular Parameterization  stretch = %f ===\n",resultCircleErr);
	if (logFile)
		fprintf(logFile,"=== Circular Parameterization  stretch = %f ===\n",resultCircleErr);

	this->resultStretch = bestStretch;
	IDtool->CleanNeighbor(BpointH,BpointT);
	return bestStretch;	
}


double    MyParameterization::SqaureParameterizationStepSampling_PARALLEL_CPU(unsigned int step_value,unsigned int &cal_count,
																			  PolarVertex *pIPV,
																			 int num_PV,
																			 FILE* logFile)
{
	double bestStretch = DBL_MAX;
	int best_startID = -1;	
	int worst_startID = -1;
	double worstStretch = -1;


	if (pU)
	{
		delete [] pU;
		pU = NULL;
	}
	if (pV)
	{
		delete [] pV;
		pV = NULL;
	}


    boundarytype=0; //square
	constsumarea3D = sqrt(1.0/sumarea3D);
	setPolarMap();
	IDList *BpointH = new IDList();
	IDList *BpointT = new IDList();
    
	BpointH->next = BpointT;
	BpointT->back = BpointH;
	
	double tlength=0.0;
	int numBorderPoint = 0;
	CalBorderPath(BpointH,BpointT,&tlength,&numBorderPoint);

	IDList *now = BpointH;
	double loop = 0;

	//ResetInnerLambda();
	setFloaterC();
	SortIndexP();
	
	//record result as array.
	vector<int> bottomLeftVertexIDList;
	
	printf("=== NUMBER BORDER EDGES : %d ===\n",numBorderPoint);
	if (logFile)
		fprintf(logFile,"=== NUMBER BORDER EDGES : %d ===\n",numBorderPoint);
	while (loop < tlength*0.25 && next(now) != BpointT)
	{
		now = next(now);		
		loop += PT->Distance(point[now->ID],point[next(now)->ID]);
		bottomLeftVertexIDList.push_back(now->ID);
	}

	//create initial A matrix  ,it is same for all condition
	//to boost up speed
	double *init_sa = NULL;
	unsigned long *init_ija = NULL;
	int nonzero=(numberV);
	for(int i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			nonzero += neighborI[i];			
		}
	}	
	init_ija = new unsigned long[2*nonzero+2];
    init_sa = new double[2*nonzero+2];	
	init_ija[1] = 2*(numberV)+2;
	int dlk=2*(numberV)+1;
  
	for(int i=0;i<numberV;i++)
	{
		init_sa[i+1] = 1.0;		
		if(boundary[i]!=1)
		{			
			PolarList *nowp = PHead[i];      
			while(nextN(nowp)!=PTail[i])
			{
				nowp = nextN(nowp);
				++dlk;
				init_sa[dlk] = -nowp->lambda;
				init_ija[dlk]= nowp->ID+1;
			}
		}
		init_ija[i+1+1]=dlk+1;
	}
	for(int i=0;i<numberV;i++)
	{
		init_sa[i+1+numberV] = 1.0;
		if(boundary[i]!=1)
		{
			PolarList *nowp = PHead[i];
			while(nextN(nowp)!=PTail[i])
			{
				nowp = nextN(nowp);
				++dlk;
				init_sa[dlk] = -nowp->lambda;
				init_ija[dlk]=nowp->ID+numberV+1;
			}
		}
		init_ija[i+numberV+1+1]=dlk+1;
	}
	//prepare for store result array
	double *bestU = NULL;
	double *bestV = NULL;
	std::vector<double *>_resultU(bottomLeftVertexIDList.size(),NULL);
	std::vector<double *>_resultV(bottomLeftVertexIDList.size(),NULL);
	std::vector<double>_resultError(bottomLeftVertexIDList.size(),0.0);
	
	int m = bottomLeftVertexIDList.size();
	cal_count = 0;
	printf("=== NUMBER of 25 %% brute force cases : %d ===\n",m);

	int step = 0;
	if (step_value  == 0 )
	{
		//calculate auto step number
		double d_step = (sqrt(numberV* ((double)m/numberF)));
		step = floor(d_step);
		if (step < 2)
			step = 2;
		printf("formulated step number: %f (%d)\n", d_step,step);
	}
	else
		step = step_value;

	
		
	int count = 1;


	#pragma omp parallel for
	for (int i = 0 ; i < bottomLeftVertexIDList.size(); i++)
	{	
		if ( ((i+1) % step) == 1 || step == 1)
		{
			#pragma omp critical
			{
				printf("S%d,",i);
			}
			SquareParametrizationCPU(bottomLeftVertexIDList[i], BpointH,BpointT, tlength,nonzero,init_sa,init_ija,_resultU[i],_resultV[i],_resultError[i],false);	
			#pragma omp critical
			{
				printf("F%d,",i);
				cal_count++;
			}
		}
		
	}
	cout << endl;

	for (int i = 0 ; i < bottomLeftVertexIDList.size(); i++)
	{
		if (_resultError[i] < bestStretch && _resultError[i] > 0)
		{
			/*
			if (bestU)
			{
				delete [] bestU;
				bestU = NULL;
			}
			if (bestV)
			{
				delete [] bestV;
				bestV= NULL;
			}
			*/

			bestU = _resultU[i];
			bestV = _resultV[i];			
			bestStretch = _resultError[i];

			best_startID = i;
		}
		else
		{			
			/*
			if (_resultU[i])
			{
				delete [] _resultU[i];
				_resultU[i]=NULL;
			}
			if (_resultV[i])
			{
				delete [] _resultV[i];
				_resultV[i] = NULL;
			}
			*/
			
		}

		if (_resultError[i] > worstStretch)
		{
			worstStretch = _resultError[i];
			worst_startID = i;
		}	
	}

	if (step > 1)
	{
		int oldbestPointIdx = best_startID;
		int neighbor_startPointIdx = best_startID;
		int neighbor_finishPointIdx = best_startID;
	
	
		for (int i = best_startID - 1; i >= 0; i --)
		{
			if (_resultError[i] <= 0)
			{
				continue;
			}
			else
			{
				neighbor_startPointIdx = i+1;
				break;
			}

		}

		for (int i = best_startID + 1; i < m; i ++)
		{
			if (_resultError[i] <= 0)
			{
				if (i == m - 1)
				{
					neighbor_finishPointIdx = i;
					break;
				}
				else
					continue;
			}
			else
			{
				neighbor_finishPointIdx = i-1;
				break;
			}

		}
		printf("=== deep check from index %d to %d\n",neighbor_startPointIdx,neighbor_finishPointIdx);
		#pragma omp parallel for
		for (int i = neighbor_startPointIdx ; i <= neighbor_finishPointIdx; i++)
		{				
			#pragma omp critical
			{
				printf("S%d,",i);
			}
			SquareParametrizationCPU(bottomLeftVertexIDList[i], BpointH,BpointT, tlength,nonzero,init_sa,init_ija,_resultU[i],_resultV[i],_resultError[i],false);	
			#pragma omp critical
			{
				printf("F%d,",i);
				cal_count++;
			}		
		}


		for (int i = 0 ; i < bottomLeftVertexIDList.size(); i++)
		{
			if (_resultError[i] < bestStretch && _resultError[i] > 0)
			{
				/*
				if (bestU)
					delete [] bestU;
				if (bestV)
					delete [] bestV;
					*/
				bestU = _resultU[i];
				bestV = _resultV[i];
			
				bestStretch = _resultError[i];
				best_startID = i;
			}
			else
			{
				/*
			
				if (_resultU[i])
				{
					delete [] _resultU[i];
					_resultU[i]=NULL;
				}
				if (_resultV[i])
				{
					delete [] _resultV[i];
					_resultV[i] = NULL;
				}
				*/
			}

			if (_resultError[i] > worstStretch)
			{
				worstStretch = _resultError[i];
				worst_startID = i;
			}	
		}


	}
	cout << endl;
	
	
	printf("=== NUMBER of 25 percent brute force cases : %d ===\n",m);
	printf("=== Find best stretch  of TEST%d  (best corner at %f) ===\n",best_startID,bestStretch);
	if (logFile)
		fprintf(logFile,"=== Find best stretch of TEST%d (best corner ERR = %f)===\n",best_startID,bestStretch);
	
	
	printf("=== Find worst stretch of TEST%d  (worst corner at %f) ===\n",worst_startID,worstStretch);
	if (logFile)
		fprintf(logFile,"=== Find worst stretch  of TEST%d  (worst corner ERR = %f) ===\n",worst_startID,worstStretch);
	

	for (int i=0;i<numberV;i++)
	{				
		pIPV[i].u = bestU[i];
		pIPV[i].v = bestV[i];
	}

	for (int i = 0 ; i < bottomLeftVertexIDList.size(); i++)
	{			
		if (_resultU[i])
		{
			delete [] _resultU[i];
			_resultU[i]=NULL;
		}
		if (_resultV[i])
		{
			delete [] _resultV[i];
			_resultV[i] = NULL;
		}
		
	}

	return bestStretch;
}


#ifdef INTEL_MKL_VERSION
double MyParameterization::ParametrizationOptimalCPU_MKL(double *ioU,double *ioV,double error,int non_zero_element,double *init_sa,unsigned long *init_ija,FILE* logFile)
{
	
	
	
	


	
#if 1
	/*---------------------------------------------------------------------------*/
  /* Define arrays for the upper triangle of the coefficient matrix and rhs vector */
  /* Compressed sparse row storage is used for sparse representation           */
  /*---------------------------------------------------------------------------*/
  MKL_INT n = 8, rci_request, itercount, expected_itercount = 8, i;
  double rhs[8];
  /* Fill all arrays containing matrix data. */
  MKL_INT ia[9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
  MKL_INT ja[18] = 
  { 1,    3,       6, 7,
       2, 3,    5,
          3,             8,
             4,       7,
                5, 6, 7,
                   6,    8,
                      7,
                         8
  };
  double a[18] = 
  { 7.E0,       1.E0,             2.E0, 7.E0,
         -4.E0, 8.E0,       2.E0,
                1.E0,                         5.E0,
                      7.E0,             9.E0,
                            5.E0, 1.E0, 5.E0,
                                 -1.E0,       5.E0,
                                       11.E0,
                                              5.E0
  };

  /*---------------------------------------------------------------------------*/
  /* Allocate storage for the solver ?par and temporary storage tmp            */
  /*---------------------------------------------------------------------------*/
  MKL_INT length = 128;
  double expected_sol[8] = { 1.E0, 0.E0, 1.E0, 0.E0, 1.E0, 0.E0, 1.E0, 0.E0 };
  /*---------------------------------------------------------------------------*/
  /* Some additional variables to use with the RCI (P)CG solver                */
  /*---------------------------------------------------------------------------*/
  double solution[8];
  MKL_INT ipar[128];
  double euclidean_norm, dpar[128], tmp[4 * 8];
  char tr = 'u';
  double eone = -1.E0;
  MKL_INT ione = 1;

  /*---------------------------------------------------------------------------*/
  /* Initialize the right hand side through matrix-vector product              */
  /*---------------------------------------------------------------------------*/
  mkl_dcsrsymv (&tr, &n, a, ia, ja, expected_sol, rhs);
  /*---------------------------------------------------------------------------*/
  /* Initialize the initial guess                                              */
  /*---------------------------------------------------------------------------*/
  for (i = 0; i < n; i++)
    solution[i] = 1.E0;
  /*---------------------------------------------------------------------------*/
  /* Initialize the solver                                                     */
  /*---------------------------------------------------------------------------*/
  dcg_init (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
  if (rci_request != 0)
    goto failure;
  /*---------------------------------------------------------------------------*/
  /* Set the desired parameters:                                               */
  /* LOGICAL parameters:                                                       */
  /* do residual stopping test                                                 */
  /* do not request for the user defined stopping test                         */
  /* DOUBLE parameters                                                         */
  /* set the relative tolerance to 1.0D-5 instead of default value 1.0D-6      */
  /*---------------------------------------------------------------------------*/
  ipar[8] = 1;
  ipar[9] = 0;
  dpar[0] = 1.E-5;
  /*---------------------------------------------------------------------------*/
  /* Check the correctness and consistency of the newly set parameters         */
  /*---------------------------------------------------------------------------*/
  dcg_check (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
  if (rci_request != 0)
    goto failure;
  /*---------------------------------------------------------------------------*/
  /* Compute the solution by RCI (P)CG solver without preconditioning          */
  /* Reverse Communications starts here                                        */
  /*---------------------------------------------------------------------------*/
rci:dcg (&n, solution, rhs, &rci_request, ipar, dpar, tmp);
  /*---------------------------------------------------------------------------*/
  /* If rci_request=0, then the solution was found with the required precision */
  /*---------------------------------------------------------------------------*/
  if (rci_request == 0)
    goto getsln;
  /*---------------------------------------------------------------------------*/
  /* If rci_request=1, then compute the vector A*tmp[0]                        */
  /* and put the result in vector tmp[n]                                       */
  /*---------------------------------------------------------------------------*/
  if (rci_request == 1)
    {
      mkl_dcsrsymv (&tr, &n, a, ia, ja, tmp, &tmp[n]);
      goto rci;
    }
  /*---------------------------------------------------------------------------*/
  /* If rci_request=anything else, then dcg subroutine failed                  */
  /* to compute the solution vector: solution[n]                               */
  /*---------------------------------------------------------------------------*/
  goto failure;
  /*---------------------------------------------------------------------------*/
  /* Reverse Communication ends here                                           */
  /* Get the current iteration number into itercount                           */
  /*---------------------------------------------------------------------------*/
getsln:dcg_get (&n, solution, rhs, &rci_request, ipar, dpar, tmp, &itercount);
  /*---------------------------------------------------------------------------*/
  /* Print solution vector: solution[n] and number of iterations: itercount    */
  /*---------------------------------------------------------------------------*/
  printf ("The system has been solved\n");
  printf ("The following solution obtained\n");
  for (i = 0; i < n / 2; i++)
    printf ("%6.3f  ", solution[i]);
  printf ("\n");
  for (i = n / 2; i < n; i++)
    printf ("%6.3f  ", solution[i]);
  printf ("\nExpected solution is\n");
  for (i = 0; i < n / 2; i++)
    {
      printf ("%6.3f  ", expected_sol[i]);
      expected_sol[i] -= solution[i];
    }
  printf ("\n");
  for (i = n / 2; i < n; i++)
    {
      printf ("%6.3f  ", expected_sol[i]);
      expected_sol[i] -= solution[i];
    }
  printf ("\nNumber of iterations: %d\n", itercount);
  i = 1;
  euclidean_norm = dnrm2 (&n, expected_sol, &i);

  /*-------------------------------------------------------------------------*/
  /* Release internal MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/
  MKL_Free_Buffers ();

  if (itercount == expected_itercount && euclidean_norm < 1.0e-12)
    {
      printf ("This example has successfully PASSED through all steps of computation!\n");
      return 0;
    }
  else
    {
      printf ("This example may have FAILED as either the number of iterations differs\n");
      printf ("from the expected number of iterations %d, or the ", expected_itercount);
      printf ("computed solution\ndiffers much from the expected solution ");
      printf ("(Euclidean norm is %e), or both.\n", euclidean_norm);
      return 1;
    }
  /*-------------------------------------------------------------------------*/
  /* Release internal MKL memory that might be used for computations         */
  /* NOTE: It is important to call the routine below to avoid memory leaks   */
  /* unless you disable MKL Memory Manager                                   */
  /*-------------------------------------------------------------------------*/
failure:printf ("This example FAILED as the solver has returned the ERROR code %d", rci_request);
  MKL_Free_Buffers ();
#endif
  return 1;
	
}
#endif

bool IsNumber(double x) 
{
    // This looks like it should always be true, 
    // but it's false if x is a NaN.
    return (x == x); 
}

double MyParameterization::ParametrizationOptimalCPU(double *ioU,double *ioV,double error,int non_zero_element,double *init_sa,unsigned long *init_ija,FILE* logFile)
{
	
	double *UaXY = new double[2*(numberV)+1];
    double *vecb = new double[2*(numberV)+1];
	
	int i;
	IDList *now = NULL;
	IDList *now2= NULL;
	PolarList *nowp=NULL;
	level=0;
  
	int nonzero=non_zero_element;


	int iter=0;
	double linerr=0.0;
	double weight=0.0;
  
	PCBCGSolver *mybcg = new PCBCGSolver(2*nonzero);
	memcpy(mybcg->sa,init_sa,sizeof(double) * (2*nonzero + 2));
	memcpy(mybcg->ija,init_ija,sizeof(unsigned long) * (2*nonzero + 2));
  

	for(i=0;i<numberV;i++)
	{    
		if(boundary[i]!=1)
		{			
			vecb[i+1] = 0.0;
			vecb[i+1+numberV] = 0.0;
		}
		else
		{			
			vecb[i+1] = ioU[i];			
			vecb[i+1+numberV] = ioV[i];
		}
		UaXY[i+1] = ioU[i];
		UaXY[i+numberV+1] = ioV[i];
	}
	int itenum = (pow((double)((numberV/20000) + 1),2)) *2000;	

	mybcg->linbcg(((unsigned long)(2*(numberV))),vecb,UaXY,1,error,itenum,&iter,&linerr);

	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			ioU[i] = UaXY[i+1];
			ioV[i] = UaXY[i+numberV+1];
		}
	}


	
	double candidate_l2 = GetStretchError(ioU,ioV);

	if (candidate_l2 < 0)
	{
		delete [] UaXY;
		delete [] vecb;
		delete mybcg;
		return candidate_l2;
	}
	double *_sigma  = new double[numberV];	
	double previous_l2= 0;
	double *sigsum = new double[numberV];
	double *prevU = new double[numberV];
	double *prevV = new double[numberV];
	
	
	do
	{
		
		previous_l2 = candidate_l2;
		memcpy(prevU,ioU,sizeof(double)*numberV);
		memcpy(prevV,ioV,sizeof(double)*numberV);
		SetSigma(ioU,ioV,_sigma);
		int dlk=2*(numberV)+1;
		for(int i=0;i<numberV;i++)
		{ 
			if(boundary[i]!=1)
			{
				sigsum[i]=0.0;
				nowp = PHead[i];
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
					++dlk;
					sigsum[i] += (-mybcg->sa[dlk]/_sigma[nowp->ID]);
				}
			}
		}
		dlk=2*(numberV)+1;
		for(int i=0;i<numberV;i++)
		{ 
			if(boundary[i]!=1)
			{
				nowp = PHead[i];
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
					++dlk;
					double newLamda = ((mybcg->sa[dlk]/_sigma[nowp->ID]))/sigsum[i];
					mybcg->sa[dlk] = newLamda;
					mybcg->sa[dlk+(nonzero - numberV)] = newLamda;

				}
			}

			UaXY[i+1] = ioU[i];
			UaXY[i+numberV+1] = ioV[i];
		}
		mybcg->linbcg(((unsigned long)(2*(numberV))),vecb,UaXY,1,error,itenum,&iter,&linerr);
		for(int i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				
				ioU[i] = UaXY[i+1];
				ioV[i] = UaXY[i+numberV+1];
					
			}
		}		
		candidate_l2 = GetStretchError(ioU,ioV);
		
 	}
	while (previous_l2 > candidate_l2 && candidate_l2 > 0);

	memcpy(ioU,prevU,sizeof(double)*numberV);
	memcpy(ioV,prevV,sizeof(double)*numberV);

	delete [] _sigma;
	delete [] sigsum;
	delete [] prevU;
	delete [] prevV;
	delete [] UaXY;
    delete [] vecb;
	delete mybcg;
	return previous_l2;

	
}


	/*
void MyParameterization::CircleParametrizationGPU(IDList *borderH,IDList *borderT,double total_length_edges, cpuCompressedMatrixType *initA, double *& resultU,double *& resultV,double &resultStretchErr,bool directsolver,bool reportlog)
{
	if (pU)
	{
		delete [] pU;
		pU = NULL;
	}
	if (pV)
	{
		delete [] pV;
		pV = NULL;
	}

	double *ioU = new double [numberV];
	double *ioV = new double [numberV];
	memset(ioU,0x00,sizeof(double)*numberV);
	memset(ioV,0x00,sizeof(double)*numberV);	

	double length = 0;
	double unitArc = (2*PI)/total_length_edges;
	IDList * now = next(borderH);
	ioU[now->ID] = 0.5*cos(0.0);
	ioV[now->ID] = 0.5*sin(0.0);
	while(next(now)!=borderT)
	{
		now = next(now);
		length += PT->Distance(point[now->ID],point[back(now)->ID]);
		//pU[now->ID] = 0.5*cos(2.0*PI*(((double)(i))/((double)(numboundary))));
		//pV[now->ID] = 0.5*sin(2.0*PI*(((double)(i))/((double)(numboundary))));
		ioU[now->ID] = 0.5*cos(unitArc*length);
		ioV[now->ID] = 0.5*sin(unitArc*length);		
	}

	
	//create initial A matrix  ,it is same for all condition	
	bool self_create_init_A = false;
	cpuCompressedMatrixType initialAMatrix(numberV*2,numberV*2);
	if (initA == NULL)
	{
		self_create_init_A = true;
		for (int i = 0; i < numberV ; i++)
		{
			if(boundary[i]!=1)
			{
			
				PolarList *nowp = PHead[i];			
				initialAMatrix(i,i) = 1.0f;
				initialAMatrix(numberV +i,numberV +i) = 1.0f;
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);						
					initialAMatrix(i,nowp->ID) = -nowp->lambda;		
					initialAMatrix(numberV +i,numberV + nowp->ID) = -nowp->lambda;
				}			
			}
			else
			{			
				//constraint
				initialAMatrix(i,i) = 1.0f;
				initialAMatrix(numberV + i,numberV +i) = 1.0f;
			}
		}
	}
	else 
		initialAMatrix = *initA;
	resultStretchErr = ParametrizationOptimalGPU(ioU,ioV,PCBCGerror,&initialAMatrix,directsolver,NULL);
	resultU =  ioU;
	resultV =  ioV;
	
	
	
}
*/
void MyParameterization::CircleParametrizationCPU(IDList *borderH,IDList *borderT,double total_length_edges,int _non_zero_element, double *_init_sa,unsigned long *_init_ija , double *& resultU,double *& resultV,double &resultStretchErr,bool reportlog)
{
	if (pU)
	{
		delete [] pU;
		pU = NULL;
	}
	if (pV)
	{
		delete [] pV;
		pV = NULL;
	}

	double *ioU = new double [numberV];
	double *ioV = new double [numberV];
	memset(ioU,0x00,sizeof(double)*numberV);
	memset(ioV,0x00,sizeof(double)*numberV);	

	double length = 0;
	double unitArc = (2*PI)/total_length_edges;
	IDList * now = next(borderH);
	ioU[now->ID] = 0.5*cos(0.0);
	ioV[now->ID] = 0.5*sin(0.0);
	while(next(now)!=borderT)
	{
		now = next(now);
		length += PT->Distance(point[now->ID],point[back(now)->ID]);
		//pU[now->ID] = 0.5*cos(2.0*PI*(((double)(i))/((double)(numboundary))));
		//pV[now->ID] = 0.5*sin(2.0*PI*(((double)(i))/((double)(numboundary))));
		ioU[now->ID] = 0.5*cos(unitArc*length);
		ioV[now->ID] = 0.5*sin(unitArc*length);		
	}

	
	double *init_sa = NULL;
	unsigned long *init_ija = NULL;
	int nonzero=(numberV);
	bool self_create_init_A = false;
	if (_init_ija == NULL || _init_sa == NULL || _non_zero_element <= 0)
	{
		self_create_init_A = true;
		for(int i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				nonzero += neighborI[i];			
			}
		}	
		init_ija = new unsigned long[2*nonzero+2];
		init_sa = new double[2*nonzero+2];	
		init_ija[1] = 2*(numberV)+2;
		int dlk=2*(numberV)+1;
  
		for(int i=0;i<numberV;i++)
		{
			init_sa[i+1] = 1.0;		
			if(boundary[i]!=1)
			{			
				PolarList *nowp = PHead[i];      
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
					++dlk;
					init_sa[dlk] = -nowp->lambda;
					init_ija[dlk]= nowp->ID+1;
				}
			}
			init_ija[i+1+1]=dlk+1;
		}
		for(int i=0;i<numberV;i++)
		{
			init_sa[i+1+numberV] = 1.0;
			if(boundary[i]!=1)
			{
				PolarList *nowp = PHead[i];
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
					++dlk;
					init_sa[dlk] = -nowp->lambda;
					init_ija[dlk]=nowp->ID+numberV+1;
				}
			}
			init_ija[i+numberV+1+1]=dlk+1;
		}
	}
	else
	{
		init_sa = _init_sa;
		init_ija = _init_ija;
		nonzero = _non_zero_element;
	}
	resultStretchErr = ParametrizationOptimalCPU(ioU,ioV,PCBCGerror,nonzero,init_sa,init_ija,NULL);
	resultU =  ioU;
	resultV =  ioV;
	
	if (self_create_init_A)
	{
		if (init_ija)
		{
			delete [] init_ija;
		}

		if (init_sa)
		{
			delete [] init_sa;
		}
	}
	
}


#if 0
double	  PsoParameterization::OptParam_Sampling(	PolarVertex *pIPV,
													int num_PV,
													FILE* logFile)
{
	printf("step: %d\n",m_step);
	double bestStretch = -1;
	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;	
	if (pU)
		delete [] pU;
	if (pV)
		delete [] pV;

	//allocate memory for u and v as number of vertex (disk topology)
	pU = new double[numberV];
	pV = new double[numberV];
	memset(pU,0x00,sizeof(double)*numberV);
	memset(pV,0x00,sizeof(double)*numberV);
	
	double *initU0 = new double[numberV];
	double *initV0 = new double[numberV];
	memset(initU0,0x00,sizeof(double)*numberV);
	memset(initV0,0x00,sizeof(double)*numberV);	

    boundarytype=0;
	setPolarMap();
	if (BpointH != NULL)
		IDtool->CleanNeighbor(BpointH,BpointT);

	BpointH = new IDList();
	BpointT = new IDList();
	

	BpointH->next = BpointT;
	BpointT->back = BpointH;
	
	total_edge_length = 0;
	int numBorderPoint = 0;
	CalBorderPath(BpointH,BpointT,&total_edge_length,&numBorderPoint);
	
	int actualParamCount = 0;
	
	totalTestCase = 0;
	actualCalCount=0;
	IDList *startPoint = BpointH->next;
	double loop = 0;
	int state=0;
	double clen=0.0;		
	
	double sum_length = 0;
	double edge_length = 0;	
	while (loop < total_edge_length*0.25)
	{		
		loop += PT->Distance(point[startPoint->ID],point[next(startPoint)->ID]);
		totalTestCase++;
		startPoint = next(startPoint);
	}
	
	if (boundary_point)
		delete [] boundary_point;
	if (boundary_point_stretch)
		delete [] boundary_point_stretch;

	boundary_point = new IDList *[totalTestCase];
	boundary_point_stretch = new double [totalTestCase];
	
	IDList *now = BpointH;
	for (int i =  0; i < totalTestCase; i++)
	{
		now = now->next;
		boundary_point[i] = now;
		boundary_point_stretch[i] = -1;
	}
	now = BpointH;	

	if (m_step  == 0 )
	{
		//calculate auto step number
		double d_step = (sqrt(numberV* ((double)totalTestCase/numberF)));
		m_step = floor(d_step);
		if (m_step < 2)
			m_step = 2;
		printf("cal auto step number: %f (%d)\n", d_step,m_step);
	}
	
		
	int count = 1;
	
	int currentBestIdx = 0;
	
	loop = 0;
	state=0;
	clen=0.0;	
	sum_length = 0;
	edge_length = 0;

	startPoint = BpointH->next;
	BpointT->ID = BpointH->next->ID; // for cal length;


	while (loop < total_edge_length*0.25)
	{
		
		state = 0;
		sum_length = 0;
		double this_side_length[4] ={0};
		int this_side_num_edge[4] = {0};		
		loop += PT->Distance(point[startPoint->ID],point[next(startPoint)->ID]);
		IDList *now = startPoint;
		if (count % m_step == 1)
		{
			do
			{				
				edge_length = PT->Distance(point[now->ID],point[next(now)->ID]);			
				if ((sum_length + edge_length>= (total_edge_length*0.25)) && state < 3)
				{
					if ((sum_length + edge_length)/(total_edge_length*0.25) >= 1.1&& this_side_num_edge[state] > 0)					
					{
						this_side_length[state] = sum_length;
						this_side_num_edge[state+1]++;
						sum_length = edge_length;
					}
					else
					{
						this_side_length[state] = sum_length+edge_length;
						this_side_num_edge[state]++;
						sum_length = 0;
					}				
					state++;
				}
				else
				{
					this_side_num_edge[state]++;
					sum_length += edge_length;
				}
				now = next(now);
				if (now == BpointT)
					now = BpointH->next;
			}
			while (now != startPoint);			
			this_side_length[state] = sum_length;


			state = 0;
			sum_length  = 0;
			clen = 0.0;
			now = startPoint;

			
			do
			{
		
				switch(state)
				{
					case 0:
						if (clen < 1.0)
						{
							pU[now->ID] = clen; 
							pV[now->ID] = 0.0;
						}
						else
						{
							pU[now->ID] = 1.0;
							pV[now->ID] = 0.0;
							sum_length = 0.0;
							state=1;											
						}
						break;
					case 1:
						if (clen < 1.0)
						{
							pU[now->ID] = 1.0; 
							pV[now->ID] = clen;
						}
						else
						{
							pU[now->ID] = 1.0;
							pV[now->ID] = 1.0;
							sum_length = 0.0;
							state=2;
						}
						break;
					case 2:
						if (clen < 1.0)
						{
							pU[now->ID] = 1.0-clen; 
							pV[now->ID] = 1.0;
						}
						else
						{
							pU[now->ID] = 0.0;
							pV[now->ID] = 1.0;
							sum_length = 0.0;
							state=3;
						}
						break;
					case 3:
						if (clen < 1.0)
						{
							pU[now->ID] = 0.0; 
							pV[now->ID] = 1.0-clen;
						}
						else
						{
							//should never enter this one;
							state=4;
						}
					default:
						break;

				}
				sum_length += PT->Distance(point[now->ID],point[next(now)->ID]);
				clen = sum_length/this_side_length[state];	
		
				now = next(now);
				if (now == BpointT)
					now = next(BpointH);
		
			}
			while(	now->ID != startPoint->ID);


		
			for(int i=0;i<numberV;i++)
			{
				if(boundary[i]==1)
				{            
		
				}
				else
				{
					pU[i] = 0;
					pV[i] = 0;
				}
			}
			ResetInnerLambda();	
		
		
			gammaP = 0.75;
			ParametrizationOptimal(iteNum,PCBCGerror,logFile);
			gammaP = 1.0;
			
			boundary_point_stretch[count - 1] = resultStretch;
			actualParamCount++;
			
			if (resultStretch < bestStretch || bestStretch < 0)
			{
				
				for (int i=0;i<numberV;i++)
				{
					pIPV[i].u = pU[i];
					pIPV[i].v = pV[i];		
				}
				

				currentBestIdx = count - 1;
				bestStretch = resultStretch;
				
			}
			
		}		
		count++;
		
		//prepare for next loop
		startPoint = next(startPoint);
		
	}



	//check neighbor of currentBest Point
	//startpoint 

	int oldbestPointIdx = currentBestIdx;
	int neighbor_startPointIdx = currentBestIdx;
	int neighbor_finishPointIdx = currentBestIdx;
	
	
	for (int i = currentBestIdx - 1; i >= 0; i --)
	{
		if (boundary_point_stretch[i] < 0)
		{
			continue;
		}
		else
		{
			neighbor_startPointIdx = i+1;
			break;
		}

	}

	for (int i = currentBestIdx + 1; i < totalTestCase; i ++)
	{
		if (boundary_point_stretch[i] < 0)
		{
			if (i == totalTestCase - 1)
			{
				neighbor_finishPointIdx = i;
				break;
			}
			else
				continue;
		}
		else
		{
			neighbor_finishPointIdx = i-1;
			break;
		}

	}


	loop = 0;
	state=0;
	clen=0.0;	
	sum_length = 0;
	edge_length = 0;

	startPoint =  boundary_point[neighbor_startPointIdx];
	int current_idx = neighbor_startPointIdx;
	
	
	while (1)
	{
		
		state = 0;
		sum_length = 0;
		double this_side_length[4] ={0};
		int this_side_num_edge[4] = {0};		
		//loop += PT->Distance(point[startPoint->ID],point[next(startPoint)->ID]);
		IDList *now = startPoint;

		if (startPoint != boundary_point[oldbestPointIdx]) //skip at oldbestPoint
		{
			do
			{				
				edge_length = PT->Distance(point[now->ID],point[next(now)->ID]);			
				if ((sum_length + edge_length>= (total_edge_length*0.25)) && state < 3)
				{
					if ((sum_length + edge_length)/(total_edge_length*0.25) >= 1.1&& this_side_num_edge[state] > 0)					
					{
						this_side_length[state] = sum_length;
						this_side_num_edge[state+1]++;
						sum_length = edge_length;
					}
					else
					{
						this_side_length[state] = sum_length+edge_length;
						this_side_num_edge[state]++;
						sum_length = 0;
					}				
					state++;
				}
				else
				{
					this_side_num_edge[state]++;
					sum_length += edge_length;
				}
				now = next(now);
				if (now == BpointT)
					now = BpointH->next;
			}
			while (now != startPoint);			
			this_side_length[state] = sum_length;


			state = 0;
			sum_length  = 0;
			clen = 0.0;
			now = startPoint;			
			do
			{
		
				switch(state)
				{
					case 0:
						if (clen < 1.0)
						{
							pU[now->ID] = clen; 
							pV[now->ID] = 0.0;
						}
						else
						{
							pU[now->ID] = 1.0;
							pV[now->ID] = 0.0;
							sum_length = 0.0;
							state=1;											
						}
						break;
					case 1:
						if (clen < 1.0)
						{
							pU[now->ID] = 1.0; 
							pV[now->ID] = clen;
						}
						else
						{
							pU[now->ID] = 1.0;
							pV[now->ID] = 1.0;
							sum_length = 0.0;
							state=2;
						}
						break;
					case 2:
						if (clen < 1.0)
						{
							pU[now->ID] = 1.0-clen; 
							pV[now->ID] = 1.0;
						}
						else
						{
							pU[now->ID] = 0.0;
							pV[now->ID] = 1.0;
							sum_length = 0.0;
							state=3;
						}
						break;
					case 3:
						if (clen < 1.0)
						{
							pU[now->ID] = 0.0; 
							pV[now->ID] = 1.0-clen;
						}
						else
						{
							//should never enter this one;
							state=4;
						}
					default:
						break;

				}
				sum_length += PT->Distance(point[now->ID],point[next(now)->ID]);
				clen = sum_length/this_side_length[state];	
		
				now = next(now);
				if (now == BpointT)
					now = next(BpointH);
		
			}
			while(	now->ID != startPoint->ID);

			for(int i=0;i<numberV;i++)
			{
				if(boundary[i]==1)
				{            
		
				}
				else
				{
					pU[i] = 0;
					pV[i] = 0;
				}
			}
			ResetInnerLambda();	
		
		
			gammaP = 0.75;
			ParametrizationOptimal(iteNum,PCBCGerror,logFile);
			gammaP = 1.0;
			
			boundary_point_stretch[current_idx] = resultStretch;
			actualParamCount++;
			
			if (resultStretch < bestStretch || bestStretch < 0)
			{
				
				for (int i=0;i<numberV;i++)
				{
					pIPV[i].u = pU[i];
					pIPV[i].v = pV[i];		
				}
				

				currentBestIdx = current_idx;
				bestStretch = resultStretch;
				
			}
			
		} //end- if (now != bestPoint)		
		current_idx++;
		
		//prepare for next loop		
		if (boundary_point[neighbor_finishPointIdx] == startPoint)
			break;  // out from loop
		else
			startPoint = next(startPoint);
	}
	



	IDtool->CleanNeighbor(BpointH,BpointT);
	
	
	printf("did (%d/%d): best is %f \n", actualParamCount,totalTestCase, bestStretch);
	return bestStretch ;
}
#endif