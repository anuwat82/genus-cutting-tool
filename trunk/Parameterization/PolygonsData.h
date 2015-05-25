#ifndef __POLYGONSDATA_HEADER__
#define __POLYGONSDATA_HEADER__

//#include "ply.h"

//#include "GL_Extension/gl_extension.h"
#include <vector>
//#include <GI/gi.h>
#include "IDCutHedge.h"
#include "MeshStructure.h"
//#include <qfileinfo.h>
#include <time.h>
//#include "gui/logwindow.h"
#include "VTK_Header.h"
#include "OpenMesh_Header.h"
#include "MyParameterization.h"
using namespace std;
#define MAX_CONNECTED_TRIANGLES_PER_VERTEX 50
#define MAX_CUT_LOOP						10

#define MODE_ORIGINAL				0x00
#define MODE_CHECKBOARD_TEXTURE		0x01
#define MODE_DOT_TEXTURE			0x02

class CPolygonsData
{
public:
	CPolygonsData(void);
	virtual ~CPolygonsData(void);
	
	void ClearData();

	//void Draw(int mode,int degree);
	//void DrawMorph(int percent);
	void	InitailDiskTopology(vtkSmartPointer<vtkPolyData> diskTopoPolydata);

	void  GetObjectCenterPosition(double *posX,double *posY,double *posZ);
	double GetObjectRadius();
	void  GetGIMfloatAdjustTransform(double *_scale,double *_translate);

	//GLuint GetGIImageTextureID();
	int Sampling(int degree,int resW, int resH, double *data_out);
	int Parameterize();
	int IteratedAugmentCut(double *calTime, vtkPolyData* op_mesh);
	int IteratedAugmentCutOriginal(double *calTime, vtkPolyData* op_mesh);
	int GetHighestCurvatureFace();
	int SquareParameterizationOptimization(unsigned int step_value, unsigned int *testCasesCount,double *calTime, vtkFloatArray *texCoord);   //step_value 1 is bruteforce   
	int SquareParameterizationManual( int mycase,double *calTime, vtkFloatArray *texCoord);  
	int CircleParameterizationOptimization(double *calTime, vtkFloatArray *texCoord, vtkDoubleArray *face_stretch); 
	int CheckBoundaryMapping(vtkDoubleArray *face_stretch);
private:
	
protected:

	//QFileInfo fi;


	//GLuint* m_uiVBO_VertexIDArray;
	//GLuint* m_uiVBO_NormalIDArray;
	//GLuint* m_uiVBO_VertexCountArray;
	unsigned int m_uiVBO_SetNumber;
	
	int most_Xindex;
	int most_minusXindex;
	int most_Yindex;
	int most_minusYindex;
	int most_Zindex;
	int most_minusZindex;

	
	
	void CalculateNormalVector();
	void FindDimension();

	//for GIs functions
	//int CreateGImeshes();
	int InitialCut();
	//int Parameterize();
	int ParameterizeOriginal();
	
	int CreateGIoutput();
	

	int  *m_connectedFace;
	void AddConnectedFace(int pointIdx,int faceIdx);	
	
	int  line_route_animation;

	double m_objectRadius;
	double m_objectCenter[3];
	
public:	
	int nelems;
	char **elist;
	int file_type;
	float version;
	int nprops;
	int num_elems;
	//PlyProperty **plist;
	Vertex *vlist;
	vector<Vertex> vertexNewCreate;
	double *vertexNormal;
	int numVertex;
	
	Face *flist;
	int *faceInfo;
	
	int			*boundarySurfaceFaceInfo[MAX_CUT_LOOP];
	int			num_boundarySurfaceFace[MAX_CUT_LOOP];
	PolarVertex *boundarySurfacePolarVertexInfo[MAX_CUT_LOOP];
	int			num_boundarySurfacePolarVertex[MAX_CUT_LOOP];


	int			*m_boundarySurfaceFaceInfo;
	int			m_num_boundarySurfaceFace;
	PolarVertex *m_boundarySurfacePolarVertexInfo;
	int			m_num_boundarySurfacePolarVertex;



	double		stretch[MAX_CUT_LOOP];
	IDCutHedge  *CutHedgeH;
	IDCutHedge  *CutHedgeT;
	int			degree_count;
	int			m_numValen2BoundaryPoint;
	bool		m_genus0closedBoundary;
	int numFace;
	
	int num_comments;
	char **comments;
	int num_obj_info;
	char **obj_info;


	int m_totalTriangles;
	int m_numBorderEdge;

	//GLs things
	//GLuint m_uiGL_giTextureID[2];
	//GLuint m_uiGL_giNormalizedTexture;

	//GIs things
	//GIuint m_uiGIMesh;
	//GIuint m_uiGIImage[3];

	static unsigned int  m_uiCheckBoardTexture;
	static int m_uiCheckBoardTextureRefCount;

	static unsigned int  m_uiDotTexture;
	static int m_uiDotTextureRefCount;


	clock_t m_calTime;
	time_t  m_startCalTime;
	time_t  m_stopCalTime;

	int SolverProcessor;  // 0 is CPU , 1 is GPU
	int SolverMethod;// 0 is iterative , 1 is direct solver
};



#endif	//__POLYGONSDATA_HEADER__
