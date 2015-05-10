#include "PolygonsData.h"

#include <math.h>
#include <memory.h>

#include "Util.h"
#include "Utils.h"
#include "stretch-minimizing.h"
#include "MyParameterization.h"
#include "customVTK\vtkFeatureEdgesEx.h"

/*
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h> 
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/Parameterization_mesh_feature_extractor.h>
typedef CGAL::Cartesian<float>     Kernel;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron3D;
typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron3D>
                                            Parameterization_polyhedron_adaptor;
*/
#define GI_RESOLUTION	(256+1)
//#define RECORD_TIME
#define FIND_BEST_CORNER

unsigned int CPolygonsData::m_uiCheckBoardTexture = 0;
int CPolygonsData::m_uiCheckBoardTextureRefCount = 0;

unsigned int CPolygonsData::m_uiDotTexture = 0;
int CPolygonsData::m_uiDotTextureRefCount = 0;
struct tempFaceStore
{
	unsigned char nverts;    /* number of vertex indices in list */
	//int *verts;              /* vertex index list */
	int verts[4];              /* vertex index list , max at 4*/	
	tempFaceStore()
	{	
		nverts = 0;
		verts[0] = verts[1] = verts[2] = verts[3] = -1;	
	}
};




CPolygonsData::CPolygonsData(void)
{
	nelems			= -1;
	elist			= NULL;

	file_type		= -1;
	version			= -1.0;

	nprops			= 0;
	
	num_elems		= 0;
	
	vlist			= NULL;
	vertexNormal	= NULL;
	numVertex		= 0;
	
	flist			= NULL;
	faceInfo		= NULL;
	numFace			= 0;
	
	
	
	
	
	num_comments	= 0;
	comments		= NULL;
	
	num_obj_info	= 0;
	obj_info		= NULL;

	m_totalTriangles = 0;




	most_Xindex			= -1;
	most_minusXindex	= -1;
	most_Yindex			= -1;
	most_minusYindex	= -1;
	most_Zindex			= -1;
	most_minusZindex	= -1;

	m_connectedFace	= NULL;

	
	
	for (int i = 0 ; i < MAX_CUT_LOOP; i++)
	{
		boundarySurfaceFaceInfo[i] = NULL;
		boundarySurfacePolarVertexInfo[i] = NULL;
		num_boundarySurfacePolarVertex[i] = 0;
		stretch[i] = -1;
	}
	m_boundarySurfaceFaceInfo = NULL;
	m_num_boundarySurfaceFace = NULL;
	m_boundarySurfacePolarVertexInfo = NULL;
	m_num_boundarySurfacePolarVertex = NULL;


	CutHedgeH = NULL;
	CutHedgeT = NULL;
	vertexNewCreate.clear();
	degree_count = 0;

	m_genus0closedBoundary = false;
	m_objectRadius  = -1;

	SolverProcessor = 0;  // 0 is CPU , 1 is GPU
	SolverMethod = 0;// 0 is iterative , 1 is direct solver
}

CPolygonsData::~CPolygonsData(void)
{
	
	ClearData();

}

void CPolygonsData::ClearData()
{		
	if (elist != NULL)
	{
		for (int i = 0 ; i < nelems; i++)
		{
			if (elist[i] != NULL)
			{
				delete [] elist[i];
				elist[i] = NULL;
			}
		}
		delete [] elist;
	}
	nelems = -1;

	file_type		= -1;
	version			= -1.0;
	nprops			= 0;
	
	num_elems		= 0;
	
	
	
	
	if (vlist != NULL)
	{		
		delete [] vlist;
	}
	if (vertexNormal != NULL)
	{
		delete [] vertexNormal;
	}
	numVertex = 0;


	if (flist != NULL)
	{	

		delete [] flist;
	}

	if (faceInfo != NULL)
		delete [] faceInfo;
	numFace = 0;

	

	if (comments != NULL)
	{
		for (int i = 0 ; i < num_comments; i++)
		{
			if (comments[i] != NULL)
			{
				delete [] comments[i];
				comments[i] = NULL;
			}
		}
		delete [] comments;
	}
	num_comments = 0;

	if (obj_info != NULL)
	{
		for (int i = 0 ; i < num_obj_info; i++)
		{
			if (obj_info[i] != NULL)
			{
				delete [] obj_info[i];
				obj_info[i] = NULL;
			}
		}
		delete [] obj_info;
	}
	num_obj_info = 0;
	m_totalTriangles = 0;


	if (m_connectedFace	!= NULL)
	{
		delete [] m_connectedFace;
		m_connectedFace = NULL;
	}


	


	

	for (int i = 0 ; i < MAX_CUT_LOOP; i++)
	{
		if (boundarySurfaceFaceInfo[i])
		{
			delete [] boundarySurfaceFaceInfo[i];
			boundarySurfaceFaceInfo[i] = NULL;
		}

		if (boundarySurfacePolarVertexInfo[i])
		{
			delete [] boundarySurfacePolarVertexInfo[i];
			boundarySurfacePolarVertexInfo[i] = NULL;
		}
	}
	
	if (CutHedgeH != NULL && CutHedgeT != NULL)
	{
		IDSet idtool ;
		idtool.CleanNeighbor(CutHedgeH,CutHedgeT);
		CutHedgeH = NULL;
		CutHedgeT = NULL;
	}
	vertexNewCreate.clear();
}








void  CPolygonsData::GetGIMfloatAdjustTransform(double *_scale,double *_translate)
{

	_translate[0] = -vlist[most_minusXindex].x;	
	_translate[1] = -vlist[most_minusYindex].y;
	_translate[2] = -vlist[most_minusZindex].z;
	

	double radius = GetObjectRadius();
	double dimX = vlist[most_Xindex].x - vlist[most_minusXindex].x;
	double dimY = vlist[most_Yindex].y - vlist[most_minusYindex].y;
	double dimZ = vlist[most_Zindex].z - vlist[most_minusZindex].z;
	
	
	double scale = 0.5f/(max(max(dimX,dimY),dimZ));
	_scale[0] = scale;
	_scale[1] = scale;
	_scale[2] = scale;
		
}

int CPolygonsData::InitialCut()
{
	
	/////**giBindMesh(m_uiGIMesh);	

	int *vertex_route = NULL;
	int *face_route = NULL;
	int num_vertex_cut_route = 0;
	clock_t calTime = clock();
	/////**giInitialCut(&vertex_route,&face_route,&num_vertex_cut_route);




	//after get route of cut-path we will split cut path to make  unboundary to be boundary surface.
	

	//to do: store information about divide edge

	IDSet idTool;
	CutHedgeH = new IDCutHedge();
	CutHedgeT = new IDCutHedge();
	CutHedgeH->next = CutHedgeT;
	CutHedgeT->back = CutHedgeH;
	IDCutHedge *pCutHedgeH =  CutHedgeH;
	IDCutHedge *pCutHedgeT =  CutHedgeT;

	//num_boundarySurfacePolarVertex[degree_count] = numVertex;
	m_num_boundarySurfacePolarVertex = numVertex;
	
	
	IDList *search = NULL; 

	

	m_numValen2BoundaryPoint = 0;
	m_numBorderEdge = num_vertex_cut_route;

	int numDuplicateVertices = 0;
	for (int i = 0 ; i < num_vertex_cut_route; i++)
	{
		
		IDCutHedge::AppendVF(vertex_route[i],face_route[i],0,pCutHedgeT);
		//Search for twin
		
		IDCutHedge *justCreated = (IDCutHedge *)pCutHedgeT->back;		
		//Search for duplicate ID		
		search = justCreated;
		while(back(search) != pCutHedgeH)
		{
			search = back(search);

			if (search->ID == justCreated->ID)
			{
				justCreated->DuplicateDegree = ((IDCutHedge *)search)->DuplicateDegree + 1;
				m_num_boundarySurfacePolarVertex++;;//num_boundarySurfacePolarVertex[degree_count]++;
				numDuplicateVertices++;
				break;
			}
		}

		if ( back(pCutHedgeT)->FaceID == back(back(pCutHedgeT))->FaceID)
		{
			m_numValen2BoundaryPoint++;
		}
	}

	

	if ( back(pCutHedgeT)->FaceID == next(pCutHedgeH)->FaceID)
	{
		m_numValen2BoundaryPoint++;
	}
	pCutHedgeT->ID = pCutHedgeH->next->ID;// last id is start point id
	
	//2013-01-08[Anuwat] store twin information
	//search->back and justCreated should be twin each other.
	IDCutHedge::RegisterTwin(pCutHedgeH,pCutHedgeT);


	if (num_vertex_cut_route == 4 && numDuplicateVertices == 1)	
		m_genus0closedBoundary = true;	
	else
		m_genus0closedBoundary = false;



	

	/*
	boundarySurfaceFaceInfo[degree_count] = new int[(numFace+(m_numValen2BoundaryPoint*2))*3 ];
	num_boundarySurfaceFace[degree_count] = numFace;

	int *p_boundarySurfaceFaceInfo = boundarySurfaceFaceInfo[degree_count];

	boundarySurfacePolarVertexInfo[degree_count] = new PolarVertex[num_boundarySurfacePolarVertex[degree_count] + m_numValen2BoundaryPoint];
	PolarVertex* p_boundarySurfacePolarVertexInfo = boundarySurfacePolarVertexInfo[degree_count];
	for (int i = 0 ; i < numVertex; i++)
	{
		p_boundarySurfacePolarVertexInfo[i].p_vertex = &vlist[i];
	}
	memcpy(p_boundarySurfaceFaceInfo,faceInfo,sizeof(int)*numFace*3);
	*/


	m_boundarySurfaceFaceInfo = new int[(numFace+(m_numValen2BoundaryPoint*2))*3 ];
	m_num_boundarySurfaceFace = numFace;
	
	m_boundarySurfacePolarVertexInfo = new PolarVertex[m_num_boundarySurfacePolarVertex + m_numValen2BoundaryPoint];	
	for (int i = 0 ; i < numVertex; i++)
	{
		m_boundarySurfacePolarVertexInfo[i].p_vertex = &vlist[i];
	}
	memcpy(m_boundarySurfaceFaceInfo,faceInfo,sizeof(int)*numFace*3);


	calTime = clock() - calTime;
	m_calTime += calTime;

	return 0;
	
}
int CPolygonsData::ParameterizeOriginal()
{
	
#ifdef RECORD_TIME
	struct tm * timeinfo = NULL;  
	timeinfo = localtime ( &m_startCalTime);
	printf ( "START CAL TIME:: %s\n", asctime (timeinfo) );
	fprintf(logFile, "START CAL TIME:: %s\n", asctime (timeinfo) );
#endif
	IDSet tool;	

	clock_t calTime = 0;
	calTime = clock();
	if ( m_genus0closedBoundary )
	{
		int *p_boundarySurfaceFaceInfo = m_boundarySurfaceFaceInfo;
		int *p_num_boundarySurfaceFaceInfo = &m_num_boundarySurfaceFace;
		PolarVertex* p_boundarySurfacePolarVertexInfo = m_boundarySurfacePolarVertexInfo;
		int *p_num_boundarySurfacePolarVertexInfo = &m_num_boundarySurfacePolarVertex;


		int num_validPolarVertex =0;		
		


		num_validPolarVertex = numVertex;
		
		MyParameterization paramTool;
		paramTool.SetPolarVertexAndFaceAndBorder(	p_boundarySurfacePolarVertexInfo,
													p_num_boundarySurfacePolarVertexInfo,
													num_validPolarVertex,
													p_boundarySurfaceFaceInfo,
													p_num_boundarySurfaceFaceInfo,
													CutHedgeH,
													CutHedgeT,
													0,
													false
												);
		int max1st = -1;
		bool boundaryFace = false;
		double circle_stretch = paramTool.GetExtremaTriangle(&max1st,&boundaryFace,1);
		int i = 0;
		int id = p_boundarySurfaceFaceInfo[max1st*3];
		while (1)
		{
			if (paramTool.neighborF[p_boundarySurfaceFaceInfo[max1st*3 + i]] >= 4)
			{
				id = p_boundarySurfaceFaceInfo[max1st*3 + i];
				break;			
			}
			else
				i++;

			if (i == 3)
			{
				max1st = paramTool.FHead[p_boundarySurfaceFaceInfo[max1st*3]]->next->ID;
				i = 0;
			}
		}
		
		
		tool.CleanNeighbor(CutHedgeH,CutHedgeT);
		CutHedgeH = new IDCutHedge();
		CutHedgeT = new IDCutHedge();
		CutHedgeH->next = CutHedgeT;
		CutHedgeT->back = CutHedgeH;

		VList *now = paramTool.VHead[id];
		while (next(now) != paramTool.VTail[id])
		{
			now = next(now);
			now = next(now);
			if (now->FaceID == max1st)
			{
				IDCutHedge::AppendVF(id,max1st,0,CutHedgeT);

				if (next(now) == paramTool.VTail[id])
					now = next(paramTool.VHead[id]);
				else
					now = next(now);

				IDCutHedge::AppendVF(now->ID,now->FaceID,0,CutHedgeT);
				now = next(now);
				int skipFace = (paramTool.neighborF[id]/2)-1;
				while (skipFace > 0)
				{
					if (next(now) == paramTool.VTail[id])
						now = paramTool.VHead[id];
					now = next(now);
					now = next(now);
					skipFace--;
				}

				IDCutHedge::AppendVF(id,now->FaceID,0,CutHedgeT);
				((IDCutHedge *)CutHedgeT->back)->DuplicateDegree = 1;
				if (next(now) == paramTool.VTail[id])
					now = next(paramTool.VHead[id]);
				else
					now = next(now);
				IDCutHedge::AppendVF(now->ID,now->FaceID,0,CutHedgeT);
				break;
			}
		}

		memcpy(p_boundarySurfaceFaceInfo,faceInfo,sizeof(int)*numFace*3);
		calTime = clock() - calTime;
		m_calTime += calTime;
	}
	// calculation
	int num_validPolarVertex =0;
	//double stopConst = 1500;
	double stopConst = 2.5;
	double previousStretch = DBL_MAX;

	int *pPrev_boundarySurfaceFaceInfo = NULL;
	int Prev_num_boundarySurfaceFaceInfo = 0;
	PolarVertex* pPrev_boundarySurfacePolarVertexInfo = NULL;
	int Prev_num_boundarySurfacePolarVertexInfo = 0;

	while (1)
	{
		bool stopLoop = false;
		calTime = clock();
		int *p_boundarySurfaceFaceInfo = m_boundarySurfaceFaceInfo;
		int *p_num_boundarySurfaceFaceInfo = &m_num_boundarySurfaceFace;
		PolarVertex* p_boundarySurfacePolarVertexInfo = m_boundarySurfacePolarVertexInfo;
		int *p_num_boundarySurfacePolarVertexInfo = &m_num_boundarySurfacePolarVertex;

		
		if (degree_count == 0)
			num_validPolarVertex = numVertex;
		
		
		MyParameterization paramTool;
		printf(	"DEGREE %d . . .\n",degree_count);
		paramTool.SetPolarVertexAndFaceAndBorder(	p_boundarySurfacePolarVertexInfo,
													p_num_boundarySurfacePolarVertexInfo,
													num_validPolarVertex,
													p_boundarySurfaceFaceInfo,
													p_num_boundarySurfaceFaceInfo,
													CutHedgeH,
													CutHedgeT,
													m_numValen2BoundaryPoint
												);	
		m_numValen2BoundaryPoint = 0;

		double current_stretch = paramTool.Parameterize(	p_boundarySurfacePolarVertexInfo,*p_num_boundarySurfacePolarVertexInfo,NULL);
		if (current_stretch <= previousStretch)
		{
			//find high stretch face
			int max1st = -1;
			bool boundaryFace = false;
			int solverMode = 0;	
			if (m_genus0closedBoundary && degree_count == 0)
			{
				solverMode = 1;
			}
			double circle_stretch = paramTool.GetExtremaTriangle(&max1st,&boundaryFace,solverMode);

			IDList *newCutPathH = NULL;
			IDList *newCutPathT = NULL;		
			int numNode = 0;			
		
			newCutPathH = new IDList();
			newCutPathT = new IDList();
			newCutPathH->next = newCutPathT;
			newCutPathT->back = newCutPathH;
			numNode = paramTool.FindShortestPathToBorderFromFace(max1st,newCutPathH,newCutPathT);
			int numDupPoint = 0;
			if (numNode < 2 )
			{			
				tool.CleanNeighbor(newCutPathH,newCutPathT);
				printf("STOP BECAUSE NUM NODE TO EXTREMA WAS %d \n",numNode);
				stopLoop = true;
			}
			else
			{			
				paramTool.AddCutBorder(CutHedgeH,CutHedgeT,newCutPathH,newCutPathT,degree_count+1,&numDupPoint);
				tool.CleanNeighbor(newCutPathH,newCutPathT);
				int new_num_boundarySurfaceFaceInfo = (*p_num_boundarySurfaceFaceInfo);		
				int new_num_boundarySurfacePolarVertexInfo = (*p_num_boundarySurfacePolarVertexInfo) + numDupPoint;

				int *new_p_boundarySurfaceFaceInfo = new int [((*p_num_boundarySurfaceFaceInfo)+4)*3];
				memcpy(new_p_boundarySurfaceFaceInfo,p_boundarySurfaceFaceInfo,sizeof(int)*(*p_num_boundarySurfaceFaceInfo)*3);

				PolarVertex* new_p_boundarySurfacePolarVertexInfo =  new PolarVertex [(*p_num_boundarySurfacePolarVertexInfo) + numDupPoint + 2];
				memcpy(new_p_boundarySurfacePolarVertexInfo,p_boundarySurfacePolarVertexInfo,sizeof(PolarVertex)*(*p_num_boundarySurfacePolarVertexInfo));

				num_validPolarVertex = (*p_num_boundarySurfacePolarVertexInfo);
				
				if (pPrev_boundarySurfaceFaceInfo)
					delete [] pPrev_boundarySurfaceFaceInfo;
				if (pPrev_boundarySurfacePolarVertexInfo)
					delete [] pPrev_boundarySurfacePolarVertexInfo;
				
				pPrev_boundarySurfaceFaceInfo = m_boundarySurfaceFaceInfo;
				Prev_num_boundarySurfaceFaceInfo = m_num_boundarySurfaceFace;
				pPrev_boundarySurfacePolarVertexInfo = m_boundarySurfacePolarVertexInfo;
				Prev_num_boundarySurfacePolarVertexInfo = m_num_boundarySurfacePolarVertex;

				m_boundarySurfaceFaceInfo = new_p_boundarySurfaceFaceInfo;
				m_num_boundarySurfaceFace = new_num_boundarySurfaceFaceInfo;
				m_boundarySurfacePolarVertexInfo = new_p_boundarySurfacePolarVertexInfo;
				m_num_boundarySurfacePolarVertex = new_num_boundarySurfacePolarVertexInfo;

				if (current_stretch >= 1.0)
					previousStretch = current_stretch;
			}
		}
		else
		{
			stopLoop = true;
		}
		
		if (stopLoop)
		{
			calTime = clock() - calTime; 
			m_calTime += calTime;
			if (pPrev_boundarySurfaceFaceInfo != NULL && pPrev_boundarySurfacePolarVertexInfo != NULL )
			{
				delete [] p_boundarySurfaceFaceInfo;
				delete [] p_boundarySurfacePolarVertexInfo;
				m_boundarySurfaceFaceInfo = pPrev_boundarySurfaceFaceInfo;
				m_num_boundarySurfaceFace = Prev_num_boundarySurfaceFaceInfo;
				m_boundarySurfacePolarVertexInfo = pPrev_boundarySurfacePolarVertexInfo;
				m_num_boundarySurfacePolarVertex = Prev_num_boundarySurfacePolarVertexInfo;
			}
			printf( "Calculate Time for seam boundary : %f sec\n", (double)m_calTime/CLOCKS_PER_SEC);
			//fprintf(logFile,"Calculate Time for seam boundary : %f sec\n", (double)m_calTime/CLOCKS_PER_SEC);			
			m_calTime = 0;
			break;			
		}

		degree_count++;
		calTime = clock() - calTime;
		m_calTime += calTime;
	}




	



	
	return 0;
}
int CPolygonsData::GetHighestCurvatureFace()
{
	if (this->numVertex == 0 || this->numFace== 0)
		return -2; //please call initial function first;

	int *p_boundarySurfaceFaceInfo = m_boundarySurfaceFaceInfo;
	int *p_num_boundarySurfaceFaceInfo = &m_num_boundarySurfaceFace;
	PolarVertex* p_boundarySurfacePolarVertexInfo = m_boundarySurfacePolarVertexInfo;
	int *p_num_boundarySurfacePolarVertexInfo = &m_num_boundarySurfacePolarVertex;

	int num_validPolarVertex =numVertex;	
	MyParameterization paramTool;
	paramTool.SetPolarVertexAndFaceAndBorder(	p_boundarySurfacePolarVertexInfo,
												p_num_boundarySurfacePolarVertexInfo,
												num_validPolarVertex,
												p_boundarySurfaceFaceInfo,
												p_num_boundarySurfaceFaceInfo,
												CutHedgeH,
												CutHedgeT,
												m_numValen2BoundaryPoint
											);	
	m_numValen2BoundaryPoint = 0;
	int max1st = -1;
	bool boundaryFace = false;
	double circle_stretch = paramTool.GetExtremaTriangle(&max1st,&boundaryFace,1);
	return max1st;
}
int CPolygonsData::Parameterize()
{	

#ifdef RECORD_TIME
	struct tm * timeinfo = NULL;  
	timeinfo = localtime ( &m_startCalTime);
	printf ( "START CAL TIME:: %s\n", asctime (timeinfo) );
	fprintf(logFile, "START CAL TIME:: %s\n", asctime (timeinfo) );
#endif
	IDSet tool;	

	clock_t calTime = 0;
	calTime = clock();
	if ( m_genus0closedBoundary )
	{
		/*		
		int *p_boundarySurfaceFaceInfo = boundarySurfaceFaceInfo[0];	
		int *p_num_boundarySurfaceFaceInfo = &num_boundarySurfaceFace[0];
		PolarVertex* p_boundarySurfacePolarVertexInfo = boundarySurfacePolarVertexInfo[0];
		int *p_num_boundarySurfacePolarVertexInfo = &num_boundarySurfacePolarVertex[0];
		*/
		int *p_boundarySurfaceFaceInfo = m_boundarySurfaceFaceInfo;
		int *p_num_boundarySurfaceFaceInfo = &m_num_boundarySurfaceFace;
		PolarVertex* p_boundarySurfacePolarVertexInfo = m_boundarySurfacePolarVertexInfo;
		int *p_num_boundarySurfacePolarVertexInfo = &m_num_boundarySurfacePolarVertex;


		int num_validPolarVertex =0;		
		


		num_validPolarVertex = numVertex;
		
		MyParameterization paramTool;
		paramTool.SetPolarVertexAndFaceAndBorder(	p_boundarySurfacePolarVertexInfo,
													p_num_boundarySurfacePolarVertexInfo,
													num_validPolarVertex,
													p_boundarySurfaceFaceInfo,
													p_num_boundarySurfaceFaceInfo,
													CutHedgeH,
													CutHedgeT,
													0,
													false
												);
		int max1st = -1;
		bool boundaryFace = false;
		double circle_stretch = paramTool.GetExtremaTriangle(&max1st,&boundaryFace,1);
		int i = 0;
		int id = p_boundarySurfaceFaceInfo[max1st*3];
		while (1)
		{
			if (paramTool.neighborF[p_boundarySurfaceFaceInfo[max1st*3 + i]] >= 4)
			{
				id = p_boundarySurfaceFaceInfo[max1st*3 + i];
				break;			
			}
			else
				i++;

			if (i == 3)
			{
				max1st = paramTool.FHead[p_boundarySurfaceFaceInfo[max1st*3]]->next->ID;
				i = 0;
			}
		}
		
		
		tool.CleanNeighbor(CutHedgeH,CutHedgeT);
		CutHedgeH = new IDCutHedge();
		CutHedgeT = new IDCutHedge();
		CutHedgeH->next = CutHedgeT;
		CutHedgeT->back = CutHedgeH;

		VList *now = paramTool.VHead[id];
		while (next(now) != paramTool.VTail[id])
		{
			now = next(now);
			now = next(now);
			if (now->FaceID == max1st)
			{
				IDCutHedge::AppendVF(id,max1st,0,CutHedgeT);

				if (next(now) == paramTool.VTail[id])
					now = next(paramTool.VHead[id]);
				else
					now = next(now);

				IDCutHedge::AppendVF(now->ID,now->FaceID,0,CutHedgeT);
				now = next(now);
				int skipFace = (paramTool.neighborF[id]/2)-1;
				while (skipFace > 0)
				{
					if (next(now) == paramTool.VTail[id])
						now = paramTool.VHead[id];
					now = next(now);
					now = next(now);
					skipFace--;
				}

				IDCutHedge::AppendVF(id,now->FaceID,0,CutHedgeT);
				((IDCutHedge *)CutHedgeT->back)->DuplicateDegree = 1;
				if (next(now) == paramTool.VTail[id])
					now = next(paramTool.VHead[id]);
				else
					now = next(now);
				IDCutHedge::AppendVF(now->ID,now->FaceID,0,CutHedgeT);
				break;
			}
		}

		IDCutHedge::RegisterTwin(CutHedgeH,CutHedgeT);
		memcpy(p_boundarySurfaceFaceInfo,faceInfo,sizeof(int)*numFace*3);
		calTime = clock() - calTime;
		m_calTime += calTime;
	}
	// calculation
	int num_validPolarVertex =0;
	//double stopConst = 1500;
	double stopConst = 2.0;
	double previousStretch = DBL_MAX;
	while (1)
	{
		bool stopLoop = false;
		calTime = clock();
		/*
		int *p_boundarySurfaceFaceInfo = boundarySurfaceFaceInfo[degree_count];	
		int *p_num_boundarySurfaceFaceInfo = &num_boundarySurfaceFace[degree_count];
		PolarVertex* p_boundarySurfacePolarVertexInfo = boundarySurfacePolarVertexInfo[degree_count];
		int *p_num_boundarySurfacePolarVertexInfo = &num_boundarySurfacePolarVertex[degree_count];
		*/
		int *p_boundarySurfaceFaceInfo = m_boundarySurfaceFaceInfo;
		int *p_num_boundarySurfaceFaceInfo = &m_num_boundarySurfaceFace;
		PolarVertex* p_boundarySurfacePolarVertexInfo = m_boundarySurfacePolarVertexInfo;
		int *p_num_boundarySurfacePolarVertexInfo = &m_num_boundarySurfacePolarVertex;

		
		if (degree_count == 0)
			num_validPolarVertex = numVertex;
		
		
		MyParameterization paramTool;
#ifdef RECORD_LOG
		printf("===== DEGREE %d =====\n",degree_count);
		fprintf(logFile,"===== DEGREE %d =====\n",degree_count);

		/*
		int numEdge = 0;
		IDCutHedge *now = CutHedgeH;
		while (next(now)!=CutHedgeT)
		{
			now = next(now);
			numEdge++;
		}
		*/
#endif
		printf(	"DEGREE %d . . .\n",degree_count);
		paramTool.SetPolarVertexAndFaceAndBorder(	p_boundarySurfacePolarVertexInfo,
													p_num_boundarySurfacePolarVertexInfo,
													num_validPolarVertex,
													p_boundarySurfaceFaceInfo,
													p_num_boundarySurfaceFaceInfo,
													CutHedgeH,
													CutHedgeT,
													m_numValen2BoundaryPoint
												);	
		m_numValen2BoundaryPoint = 0;
#if 0
		if (degree_count >= 10)
		{
			char str[256] = "";
			sprintf(str,"%s_degree%d.ply2",fi.fileName().toAscii().data(),degree_count);
			FILE *out = fopen(str,"wt");
			fprintf(out,"%d\n"
							"%d\n",*p_num_boundarySurfacePolarVertexInfo,*p_num_boundarySurfaceFaceInfo);
			for (int i = 0 ; i < *p_num_boundarySurfacePolarVertexInfo; i++)
			{
				fprintf(out,"%f\n"
							"%f\n"
							"%f\n", p_boundarySurfacePolarVertexInfo[i].p_vertex->x,
									p_boundarySurfacePolarVertexInfo[i].p_vertex->y,
									p_boundarySurfacePolarVertexInfo[i].p_vertex->z);
			}
			for (int i = 0 ; i < *p_num_boundarySurfaceFaceInfo; i++)
			{									
				fprintf(out,"3\n"
								"%d\n"
								"%d\n"
								"%d\n",	p_boundarySurfaceFaceInfo[i*3 + 0],											
										p_boundarySurfaceFaceInfo[i*3 + 1],
										p_boundarySurfaceFaceInfo[i*3 + 2]);
			}
			fclose(out);
		}		
#endif	
#ifdef RECORD_LOG
		printf("NUMBER POLAR VERTEX: %d\n",*p_num_boundarySurfacePolarVertexInfo);
		fprintf(logFile,"NUMBER POLAR VERTEX: %d\n",*p_num_boundarySurfacePolarVertexInfo);
		printf("NUMBER SURFACE TRIANGLES: %d\n",*p_num_boundarySurfaceFaceInfo);
		fprintf(logFile,"NUMBER SURFACE TRIANGLES: %d\n",*p_num_boundarySurfaceFaceInfo);
		
		double idealEdgeLength = sqrt((*p_num_boundarySurfaceFaceInfo)*8.0);
		printf("NUMBER OF EDGE BOUNDARY: %d (%d)\n",numEdge,(int)idealEdgeLength);
		fprintf(logFile,"NUMBER OF EDGE BOUNDARY: %d (%d)\n",numEdge,(int)idealEdgeLength);		
#endif	
		
#ifdef RECORD_LOG				
		printf("SQUARE STRETCH: %f\n",stretch[degree_count]);
		fprintf(logFile,"SQUARE STRETCH: %f\n",stretch[degree_count]);
		//break;
#endif
		
		int max1st = -1;
		bool boundaryFace = false;
		int solverMode = 0;
		int numDupPoint = 0;

		if (m_genus0closedBoundary && degree_count == 0)
		{
			solverMode = 1;
		}
		double circle_stretch = 0;
		circle_stretch = paramTool.GetExtremaTriangle(&max1st,&boundaryFace,solverMode);
#if 0

		fprintf(logFile,"FLOATER CIRCLE STRETCH: %f\n",circle_stretch);
		fprintf(logFile,"MAX INNER FACE STRETCH: %d\n",max1st);
		fprintf(logFile,"\n");

		printf(	"FLOATER CIRCLE STRETCH: %f\n"
				"MAX INNER FACE STRETCH: %d\n"
				"\n",circle_stretch,max1st);
#endif		
		//break;
		/*
		if (boundaryFace && degree_count == 13)
		{
			boundaryFace = false;
			max1st = 66549;
		}
		*/
		
		if (boundaryFace )
		{
			printf(	"STOP BECAUSE BOUNDARY FACE (%f/%f)\n",circle_stretch,stopConst);
			//fprintf(logFile,	"STOP BECAUSE BOUNDARY FACE (%f/%f)\n",circle_stretch,stopConst);
			calTime = clock() - calTime; 
			m_calTime += calTime;
			stopLoop = true;

			
		}
		else if (circle_stretch < stopConst)
		{
			printf(	"STOP BECAUSE LOWER THAN CONST (%f/%f)\n",circle_stretch,stopConst);
			//fprintf(logFile,"STOP BECAUSE LOWER THAN CONST (%f/%f)\n",circle_stretch,stopConst);
			calTime = clock() - calTime; 
			m_calTime += calTime;
			stopLoop = true;
		}
		else
		{
			IDList *newCutPathH = NULL;
			IDList *newCutPathT = NULL;		
			int numNode = 0;			
		
			newCutPathH = new IDList();
			newCutPathT = new IDList();
			newCutPathH->next = newCutPathT;
			newCutPathT->back = newCutPathH;
			numNode = paramTool.FindShortestPathToBorderFromFace(max1st,newCutPathH,newCutPathT);		
		
		
			if (numNode < 4 )
			{			
				tool.CleanNeighbor(newCutPathH,newCutPathT);
				printf("STOP BECAUSE NUM NODE TO EXTREMA WAS %d (%f/%f)\n",numNode,circle_stretch,stopConst);
				//fprintf(logFile,"STOP BECAUSE NUM NODE TO EXTREMA WAS %d (%f/%f)\n",numNode,circle_stretch,stopConst);
				calTime = clock() - calTime; 
				m_calTime += calTime;
				stopLoop = true;
			}
			else
			{			
				paramTool.AddCutBorder(CutHedgeH,CutHedgeT,newCutPathH,newCutPathT,degree_count+1,&numDupPoint);
				tool.CleanNeighbor(newCutPathH,newCutPathT);
			
			}
		}
		stopLoop = true;
		if (stopLoop)
		{
#ifdef RECORD_TIME
			time(&m_stopCalTime);
			timeinfo = localtime ( &m_stopCalTime);
			printf ( "STOP CAL TIME:: %s\n", asctime (timeinfo) );
			fprintf(logFile, "STOP CAL TIME:: %s\n", asctime (timeinfo) );
#endif
			
			printf( "Calculate Time for seam boundary : %f sec\n", (double)m_calTime/CLOCKS_PER_SEC);
			//fprintf(logFile,"Calculate Time for seam boundary : %f sec\n", (double)m_calTime/CLOCKS_PER_SEC);			
			m_calTime = 0;
			
#ifdef FIND_BEST_CORNER
#ifdef RECORD_TIME			
			time(&m_startCalTime);
			timeinfo = localtime ( &m_startCalTime);
			printf ( "START CAL TIME:: %s\n", asctime (timeinfo) );
			fprintf(logFile, "START CAL TIME:: %s\n", asctime (timeinfo) );
#endif
			calTime = clock();
			//paramTool.PARAM_MYEXPER4(	p_boundarySurfacePolarVertexInfo,*p_num_boundarySurfacePolarVertexInfo,logFile);
			//paramTool.Parameterize(	p_boundarySurfacePolarVertexInfo,*p_num_boundarySurfacePolarVertexInfo,logFile);
			//paramTool.PARAM_PARALLEL_GPU(	p_boundarySurfacePolarVertexInfo,*p_num_boundarySurfacePolarVertexInfo,logFile);
			paramTool.PARAM_PARALLEL_CPU(	p_boundarySurfacePolarVertexInfo,*p_num_boundarySurfacePolarVertexInfo,NULL);
			calTime = clock() - calTime; 
			m_calTime += calTime;
			//printf( "Calculate Time for seam boundary : %f sec\n", (double)m_calTime/CLOCKS_PER_SEC);
			//fprintf(logFile,"Calculate Time for seam boundary : %f sec\n", (double)m_calTime/CLOCKS_PER_SEC);			
			
			printf( "Calculate Time for corner : %f sec\n", (double)m_calTime/CLOCKS_PER_SEC);
			//fprintf(N,"Calculate Time for corner : %f sec\n", (double)m_calTime/CLOCKS_PER_SEC);


#ifdef RECORD_TIME
			time(&m_stopCalTime);
			timeinfo = localtime ( &m_stopCalTime);
			printf ( "STOP CAL TIME:: %s\n", asctime (timeinfo) );
			fprintf(logFile, "STOP CAL TIME:: %s\n", asctime (timeinfo) );
#endif
#endif
			break;
		}
		
		int new_num_boundarySurfaceFaceInfo = (*p_num_boundarySurfaceFaceInfo);		
		int new_num_boundarySurfacePolarVertexInfo = (*p_num_boundarySurfacePolarVertexInfo) + numDupPoint;

		int *new_p_boundarySurfaceFaceInfo = new int [((*p_num_boundarySurfaceFaceInfo)+4)*3];
		memcpy(new_p_boundarySurfaceFaceInfo,p_boundarySurfaceFaceInfo,sizeof(int)*(*p_num_boundarySurfaceFaceInfo)*3);

		PolarVertex* new_p_boundarySurfacePolarVertexInfo =  new PolarVertex [(*p_num_boundarySurfacePolarVertexInfo) + numDupPoint + 2];
		memcpy(new_p_boundarySurfacePolarVertexInfo,p_boundarySurfacePolarVertexInfo,sizeof(PolarVertex)*(*p_num_boundarySurfacePolarVertexInfo));

		/*
		num_boundarySurfacePolarVertex[degree_count+ 1] =  (*p_num_boundarySurfacePolarVertexInfo) + numDupPoint;
		num_boundarySurfaceFace[degree_count+1] = (*p_num_boundarySurfaceFaceInfo);

		boundarySurfaceFaceInfo[degree_count+1] = new int [((*p_num_boundarySurfaceFaceInfo)+4)*3];
		memcpy(boundarySurfaceFaceInfo[degree_count+1],p_boundarySurfaceFaceInfo,sizeof(int)*(*p_num_boundarySurfaceFaceInfo)*3);
		
		boundarySurfacePolarVertexInfo[degree_count+1] = new PolarVertex [(*p_num_boundarySurfacePolarVertexInfo) + numDupPoint + 2];
		memcpy(boundarySurfacePolarVertexInfo[degree_count+1],p_boundarySurfacePolarVertexInfo,sizeof(PolarVertex)*(*p_num_boundarySurfacePolarVertexInfo));
		*/

		num_validPolarVertex = (*p_num_boundarySurfacePolarVertexInfo);
		delete [] p_boundarySurfaceFaceInfo;
		delete [] p_boundarySurfacePolarVertexInfo;
		m_boundarySurfaceFaceInfo = new_p_boundarySurfaceFaceInfo;
		m_num_boundarySurfaceFace = new_num_boundarySurfaceFaceInfo;
		m_boundarySurfacePolarVertexInfo = new_p_boundarySurfacePolarVertexInfo;
		m_num_boundarySurfacePolarVertex = new_num_boundarySurfacePolarVertexInfo;


		//previousStretch = stretch[degree_count];
		degree_count++;
		calTime = clock() - calTime;
		m_calTime += calTime;
	}

	

	//IDCutHedge::ReportTwin(CutHedgeH,CutHedgeT,m_boundarySurfacePolarVertexInfo);


	return 0;
	
}



int CPolygonsData::Sampling(int degree,int resW, int resH, double *data_out)
{
	/*
	if (degree > degree_count)
	{
		return -1;
	}
	*/
	if (resW <= 0 || resH <= 0)
	{
		return -1;
	}
	//PolarVertex *polar	=  boundarySurfacePolarVertexInfo[degree];
	//int			*face	= boundarySurfaceFaceInfo[degree];

	PolarVertex *polar	=  m_boundarySurfacePolarVertexInfo;
	int			*face	= m_boundarySurfaceFaceInfo;

	double Width  = resW;
	double Height = resH;
	double sampling_u;
	double sampling_v;
	for (int i = 0 ; i < resH; i++)
	{
		sampling_v = (double)i/(double)(resH-1);
		double T = (double )i;
		for (int j = 0 ; j < resW; j++)
		{
			//get triangle that  this sampling point is inside.
			double lamda0,lamda1,lamda2;
			sampling_u = (double)j/(double)(resW-1);
			int k = 0;
			for (k = 0; k < m_num_boundarySurfaceFace; k++)
			{
				
				if(( min(min(polar[face[k*3 + 0]].u,polar[face[k*3 + 1]].u),polar[face[k*3 + 2]].u) <= sampling_u)&&			
				   ( max(max(polar[face[k*3 + 0]].u,polar[face[k*3 + 1]].u),polar[face[k*3 + 2]].u) >= sampling_u)&&
				   ( min(min(polar[face[k*3 + 0]].v,polar[face[k*3 + 1]].v),polar[face[k*3 + 2]].v) <= sampling_v)&&			
				   ( max(max(polar[face[k*3 + 0]].v,polar[face[k*3 + 1]].v),polar[face[k*3 + 2]].v) >= sampling_v))
				
				{
					if (isPointInsideTriangle(	
									sampling_u,sampling_v,
									polar[face[k*3 + 0]].u,polar[face[k*3 + 0]].v,
									polar[face[k*3 + 1]].u,polar[face[k*3 + 1]].v,
									polar[face[k*3 + 2]].u,polar[face[k*3 + 2]].v,
									&lamda0,
									&lamda1,
									&lamda2))
						break;
				}				
			}
			
			if (k == m_num_boundarySurfaceFace)
				return -1;			
			data_out[(resW*3*i) + j*3 + 0] =	lamda0*polar[face[k*3 + 0]].p_vertex->x + 
												lamda1*polar[face[k*3 + 1]].p_vertex->x + 
												lamda2*polar[face[k*3 + 2]].p_vertex->x;
			
			data_out[(resW*3*i) + j*3 + 1] =	lamda0*polar[face[k*3 + 0]].p_vertex->y + 
												lamda1*polar[face[k*3 + 1]].p_vertex->y + 
												lamda2*polar[face[k*3 + 2]].p_vertex->y;

			data_out[(resW*3*i) + j*3 + 2] =	lamda0*polar[face[k*3 + 0]].p_vertex->z + 
												lamda1*polar[face[k*3 + 1]].p_vertex->z + 
												lamda2*polar[face[k*3 + 2]].p_vertex->z;
			
		}
	}

	return 0;
}


int CPolygonsData::CreateGIoutput()
{
	/*
	int res = GI_RESOLUTION;

	
	int iAttribCount, iNRes = (res&1) ? (res<<1)-1 : (res<<1);

	
	glGenTextures(2, m_uiGL_giTextureID);
	glBindTexture(GL_TEXTURE_2D, m_uiGL_giTextureID[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, res, res, 0, GL_RGBA, GL_FLOAT, NULL);
	glBindTexture(GL_TEXTURE_2D, m_uiGL_giTextureID[1]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, iNRes, iNRes, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

	// create images and sample mesh attributes into them
	giGenImages(3, m_uiGIImage);
	giBindImage(GI_GEOMETRY,	m_uiGIImage[0]);
	giBindImage(GI_NORMAL,		m_uiGIImage[1]);
	giBindImage(GI_TEXTURE,		m_uiGIImage[2]);

	
	if (m_uiCheckBoardTexture == 0 && m_uiCheckBoardTextureRefCount == 0)
	{
		unsigned char checkerboard[256][256];
		// make texture for param visiualization
		for (int i=0; i<256*4 + 1; ++i)
			for(int j=0; j<256*4 + 1; ++j)
		{
			if (i%4 == 0 && j%4 == 0)
				checkerboard[i][j] = 255;
			else
				checkerboard[i][j] = 0;
		}
		glGenTextures(1, &m_uiCheckBoardTexture);
		glBindTexture(GL_TEXTURE_2D, m_uiCheckBoardTexture);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE,256 , 256, 0, 
			GL_LUMINANCE, GL_UNSIGNED_BYTE, checkerboard);
	}
	m_uiCheckBoardTextureRefCount++;

	giImageGLTextureData(GI_GEOMETRY, res, res, 4, GL_FLOAT, m_uiGL_giTextureID[0]);
	giImageGLTextureData(GI_NORMAL, iNRes, iNRes, 3, GI_UNSIGNED_BYTE, m_uiGL_giTextureID[1]);
	giImageGLTextureData(GI_TEXTURE, 256, 256, 1, GI_UNSIGNED_BYTE, m_uiCheckBoardTexture);
	giSamplerParameteri(GI_SAMPLER, GI_SAMPLER_SOFTWARE);
	giSample(GI_GEOMETRY_BIT | GI_NORMAL_BIT);
	giGetIntegerv(GI_SAMPLED_ATTRIB_COUNT, &iAttribCount);
	if(iAttribCount != 2)
	{
		giDeleteImages(3, m_uiGIImage);
		return -1;
	}
	else
	{

		
	}
	*/
	return 0;
}

void CPolygonsData::GetObjectCenterPosition(double *posX,double *posY,double *posZ)
{
	if (numVertex > 0)
	{
		*posX = (vlist[most_Xindex].x + vlist[most_minusXindex].x)/2.0;
		*posY = (vlist[most_Yindex].y + vlist[most_minusYindex].y)/2.0;
		*posZ = (vlist[most_Zindex].z + vlist[most_minusZindex].z)/2.0;

		
	}
	else
	{
		*posX = 0.0;
		*posY = 0.0;
		*posZ = 0.0;
	}

}


double CPolygonsData::GetObjectRadius()
{
	if (numVertex > 0)
	{
		Vertex center;
		GetObjectCenterPosition(&center.x,&center.y,&center.z);	
		
		double v[3] = {	vlist[most_Xindex].x - center.x,
						vlist[most_Yindex].y - center.y,
						vlist[most_Zindex].z - center.z};
		double length = sqrt(
								pow(v[0],2) +
								pow(v[1],2) +
								pow(v[2],2)
							);
		m_objectRadius = length;
		return length;
	}
	else
	{
		return 0.0;
	}
}


void CPolygonsData::FindDimension()
{
	if (numVertex > 0)
	{
		most_Xindex			= 0;
		most_minusXindex	= 0;
		most_Yindex			= 0;
		most_minusYindex	= 0;
		most_Zindex			= 0;
		most_minusZindex	= 0;
		
		for (int i = 1 ; i < numVertex ; i++)	//start from 1
		{
			//X
			if (vlist[i].x > vlist[most_Xindex].x)
			{
				most_Xindex = i;
			}
			else if (vlist[i].x < vlist[most_minusXindex].x)
			{
				most_minusXindex = i;
			}

			//Y
			if (vlist[i].y > vlist[most_Yindex].y)
			{
				most_Yindex = i;
			}
			else if (vlist[i].y < vlist[most_minusYindex].y)
			{
				most_minusYindex = i;
			}

			//Z
			if (vlist[i].z > vlist[most_Zindex].z)
			{
				most_Zindex = i;
			}
			else if (vlist[i].z < vlist[most_minusZindex].z)
			{
				most_minusZindex = i;
			}
		}
	}
}

void CPolygonsData::CalculateNormalVector()
{
	if (vertexNormal!= NULL)
	{
		delete [] vertexNormal;
		vertexNormal = NULL;
	}
	vertexNormal	= new double[numVertex*3];

	for (int face_idx = 0 ; face_idx < numFace; face_idx++)
	{
		switch(flist[face_idx].nverts)
		{
			case 3:
				{
					int vertex_idx[3] = {
											flist[face_idx].faceInfo_ptr[0],
											flist[face_idx].faceInfo_ptr[1],
											flist[face_idx].faceInfo_ptr[2],
										};
					double v[3] = {	vlist[vertex_idx[1]].x - vlist[vertex_idx[0]].x,
									vlist[vertex_idx[1]].y	- vlist[vertex_idx[0]].y,
									vlist[vertex_idx[1]].z	- vlist[vertex_idx[0]].z };
					double w[3] = {	vlist[vertex_idx[2]].x - vlist[vertex_idx[0]].x,
									vlist[vertex_idx[2]].y	- vlist[vertex_idx[0]].y,
									vlist[vertex_idx[2]].z	- vlist[vertex_idx[0]].z };


					double CrossProduct[3] = {	v[1]*w[2] - w[1]*v[2],
												w[0]*v[2] - v[0]*w[2],
												v[0]*w[1] - w[0]*v[1] };

					vector_normalize(CrossProduct, CrossProduct);
					vector_copy(CrossProduct,flist[face_idx].normalVector);
				}
				break;
			case 1:
				{
					int vertex_idx = flist[face_idx].faceInfo_ptr[0];
					double v[3] = {	vlist[vertex_idx].x,
									vlist[vertex_idx].y,
									vlist[vertex_idx].z };
					vector_normalize(v, &vertexNormal[vertex_idx*3]);
				}
				break;	

		}
		
	}

	if (m_connectedFace != NULL)
	{

		for (int vertex_idx = 0 ; vertex_idx < numVertex; vertex_idx++)
		{
			double v[3] = {0};
			double n[3] = {0};

			for (int i = 0 ; i < MAX_CONNECTED_TRIANGLES_PER_VERTEX  ; i++)
			{
				int fi = m_connectedFace[(vertex_idx*MAX_CONNECTED_TRIANGLES_PER_VERTEX) + i];
				if (fi < 0)
					break;
				
				vector_copy(flist[fi].normalVector,n);
				if (vector_length(n) > 0)
				{
					if (vector_length(v) > 0)
					{					
						double dotResult = vector_dot(v,n);
						if (!(dotResult > 0.707)) // not less than 45 degree
						{
							vector_scalarMultiply(n,-1.0f,n);
							dotResult = vector_dot(v,n);
							if (!(dotResult > 0.707)) // not less than 45 degree
							{							
								continue;
							}
						}
					}
					vector_add(v,n,v);
					vector_normalize(v,v);
				}
			}
			vector_copy(v,&vertexNormal[vertex_idx*3]);
		}
	}
	
}









void CPolygonsData::AddConnectedFace(int pointIdx,int faceIdx)
{
	if (m_connectedFace == NULL)
	{
		return;
	}
	for(int i = 0 ; i < MAX_CONNECTED_TRIANGLES_PER_VERTEX ; i++)
	{
		if (m_connectedFace[(pointIdx*MAX_CONNECTED_TRIANGLES_PER_VERTEX)+i] < 0)
		{
			m_connectedFace[(pointIdx*MAX_CONNECTED_TRIANGLES_PER_VERTEX)+i] = faceIdx;
			return;
		}
	}
	printf("Error: MAX_CONNECTED_TRIANGLES_PER_VERTEX exceeded from vertex#%d\n",pointIdx);
}

void	CPolygonsData::InitailDiskTopology(vtkSmartPointer<vtkPolyData> diskTopoPolydata)
{
	OmMesh mesh;
	vtkPolydata2OpenMesh(diskTopoPolydata,&mesh );	
	int _numV =  mesh.n_vertices();
	int bid ;
	for ( bid = _numV -1 ; bid >= 0; bid--)
	{
		if (mesh.is_boundary (mesh.vertex_handle(bid)))
			break;
	}

	//store vlist and flist
	int numV = numVertex = diskTopoPolydata->GetNumberOfPoints();
	int numF = numFace = diskTopoPolydata->GetNumberOfPolys();
	
	vlist = (Vertex *) new Vertex[numV];
	for (int i = 0 ; i < numV ; i++)
	{
		vlist[i].Set(diskTopoPolydata->GetPoint(i));
	}
	m_num_boundarySurfacePolarVertex = numVertex;
	m_numValen2BoundaryPoint = 0;
	m_numBorderEdge = 0;
	OmMesh::HalfedgeHandle firstBorderHE = mesh.halfedge_handle( mesh.vertex_handle(bid));
	OmMesh::HalfedgeHandle currentBorderHE =firstBorderHE;
	IDSet idTool;
	CutHedgeH = new IDCutHedge();
	CutHedgeT = new IDCutHedge();
	CutHedgeH->next = CutHedgeT;
	CutHedgeT->back = CutHedgeH;
	IDCutHedge *pCutHedgeH =  CutHedgeH;
	IDCutHedge *pCutHedgeT =  CutHedgeT;

	do
	{
		int vid = mesh.from_vertex_handle(currentBorderHE).idx();		
		int fid = mesh.opposite_face_handle(currentBorderHE).idx();
		IDCutHedge::AppendVF(vid,fid,0,pCutHedgeT);
		m_numBorderEdge++;
		if ( back(pCutHedgeT)->FaceID == back(back(pCutHedgeT))->FaceID)
		{
			m_numValen2BoundaryPoint++;
		}
		currentBorderHE = mesh.next_halfedge_handle(currentBorderHE);
	}
	while (currentBorderHE!= firstBorderHE);
	
	if ( back(pCutHedgeT)->FaceID == next(pCutHedgeH)->FaceID)
	{
		m_numValen2BoundaryPoint++;
	}
	pCutHedgeT->ID = pCutHedgeH->next->ID;// last id is start point id


	m_boundarySurfaceFaceInfo = new int[(numFace+(m_numValen2BoundaryPoint*2))*3 ];
	m_num_boundarySurfaceFace = numFace;
	
	m_boundarySurfacePolarVertexInfo = new PolarVertex[m_num_boundarySurfacePolarVertex + m_numValen2BoundaryPoint];	
	for (int i = 0 ; i < numVertex; i++)
	{
		m_boundarySurfacePolarVertexInfo[i].p_vertex = &vlist[i];
	}

	diskTopoPolydata->GetPolys()->InitTraversal();
	vtkIdType npts;
	vtkIdType *pointID;
	int *storeFace = m_boundarySurfaceFaceInfo;
	while(diskTopoPolydata->GetPolys()->GetNextCell(npts,pointID) != 0)
	{
		if (npts != 3)
			throw;
		storeFace[0] = pointID[0];
		storeFace[1] = pointID[1];
		storeFace[2] = pointID[2];
		storeFace += 3;
	}


	return;
}