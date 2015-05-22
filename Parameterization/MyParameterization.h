#pragma once

#include <math.h>
#include "Polyhedron.h"
#include "IDCutHedge.h"
#include "PolygonsData.h"


#include <mkl.h>    //comment out this line if do not have intel mkl
#ifdef INTEL_MKL_VERSION 
#include <boost/numeric/ublas/matrix_sparse.hpp>
typedef boost::numeric::ublas::compressed_matrix<double,boost::numeric::ublas::row_major,0,boost::numeric::ublas::unbounded_array<MKL_INT>>       cpuCompressedMatrixTypeIndex0;
typedef boost::numeric::ublas::compressed_matrix<double,boost::numeric::ublas::row_major,1,boost::numeric::ublas::unbounded_array<MKL_INT>>       cpuCompressedMatrixTypeIndex1;
#else
#include <boost/numeric/ublas/matrix_sparse.hpp>
typedef boost::numeric::ublas::compressed_matrix<double,boost::numeric::ublas::row_major,0,boost::numeric::ublas::unbounded_array<int>>       cpuCompressedMatrixTypeIndex0;
typedef boost::numeric::ublas::compressed_matrix<double,boost::numeric::ublas::row_major,1,boost::numeric::ublas::unbounded_array<int>>       cpuCompressedMatrixTypeIndex1;
#endif
typedef struct MirrorFace 
{
	Point3d pt[3];	
	double area3d;
	MirrorFace()
	{		
		area3d = 0;
	}
} MirrorFace;

class MyParameterization : public Polyhedron
{
public:
	MyParameterization(void);
	virtual ~MyParameterization(void);
	
	


	void    SetPolarVertexAndFaceAndBorder(	PolarVertex *pIPV,
											int *pIOnum_PV,
											int num_originalVertex,
											int	*pIOFace,
											int *pIOnum_Face,
											IDCutHedge *pCutHedgeH,
											IDCutHedge *pCutHedgeT,
											int num_valen2Boundary,
											bool handleValence2 = true);
	void	AddCutBorder(						
							IDCutHedge *cutHedgeH,
							IDCutHedge *cutHedgeT,
							IDList *pAddHead,
							IDList *pAddTail,
							int degree,
							int *numDuplicate
						);

	double	Parameterize(PolarVertex *pIPV,
						 int num_PV,
						 FILE* logFile= NULL);
	double	CircularParameterize(PolarVertex *pIPV,
								 int num_PV,
								 FILE* logFile= NULL);
	double	CircularParameterizeOptimalEx(PolarVertex *pIPV,
										 int num_PV,
										 FILE* logFile= NULL);

	void	linbcg_Solve();
	#ifdef INTEL_MKL_VERSION
	void	mkl_Solve();
	#endif
	void	DeleteFace(int id);
	void	ChangeVertexIndexInSurfaceFaces(int oldID, int newID, int faceExcludeID, int startHalfEdge,int stopHalfEdge,int	*pIOFace = NULL); 
	void	CalculateEdgeLength();
	int		FindShortestPathToBorderFromPoint( int startVertexID,IDList *pHead,IDList *pTail,double *length = NULL);
	int		FindShortestPathToBorderFromFace( int startFaceID,IDList *pHead,IDList *pTail);
	bool    InsertDijkstra(int ID, double length,IDList *pHead,IDList *pTail);
	double	GetExtremaTriangle(int *index , bool *borderFace ,int SolverMode = 0);
	double	GetExtremaVertex(int *index , bool *borderVertex, int SolverMode = 0);
	void    setTutteC();
	void	MyBoundaryMap();
	bool	HandleValence2Boundary(
									PolarVertex *pIPV,
									int *pIOnum_PV,
									int	*pIOFace,
									int *pIOnum_Face,									
									IDCutHedge *pCutHedgeH,
									IDCutHedge *pCutHedgeT
								  );
	void	SortV();
	double	ParamWithFaceNormalStretch(PolarVertex *pIPV,
									 int num_PV,
									 FILE* logFile= NULL);
	
	double	DAM_Param(	PolarVertex *pIPV,
						int num_PV,
						FILE* logFile= NULL);
	void	setDAMC();
	void	setSigmaArea();

	void	createMirrorFaceData();
	double    PARAM_MYEXPER(PolarVertex *pIPV,
									 int num_PV,
									 FILE* logFile= NULL);
	double    PARAM_MYEXPER2(PolarVertex *pIPV,
									 int num_PV,
									 FILE* logFile= NULL);

	double    PARAM_MYEXPER3(PolarVertex *pIPV,
							 int num_PV,
							 FILE* logFile= NULL);

	double    PARAM_MYEXPER4(PolarVertex *pIPV,
							 int num_PV,
							 FILE* logFile= NULL);
	/*
	double    PARAM_VIENNA_CL(PolarVertex *pIPV,
							 int num_PV,
							 FILE* logFile= NULL);
	double    PARAM_PARALLEL_GPU(PolarVertex *pIPV,
							 int num_PV,
							 FILE* logFile= NULL);
*/
	double    PARAM_PARALLEL_CPU(PolarVertex *pIPV,
							 int num_PV,
							 FILE* logFile= NULL);


	double    SqaureParameterizationStepSampling_PARALLEL_CPU( unsigned int step_value,
															   unsigned int &cal_count,
																PolarVertex *pIPV,
															 int num_PV,
															 FILE* logFile= NULL);

	void	ResetInnerLambda();
	void    CalBorderPath(IDList *BPointH,IDList *BPointT,double *length,int *numPoint);

	double	getCurrentE_EX();
	double	getCurrentE_U();
	double	getCurrentE_V();
	void	setPolarMap_EX();


	void	setSigmaU();
	void	setSigmaV();
	

	void    FindMinMaxStretchBoundaryFace(int *minID,int *maxID);


	void ParametrizationSmoothOptimal_EX(double startGamma, double finishGamma,int itenum,double error,FILE* logFile = NULL);
	void ParametrizationOptimalSaveU0(double *op_u0, double *op_v0,int itenum,double error,FILE* logFile = NULL);
	//void ParametrizationOptimalGPU(double error,FILE* logFile);
	//double ParametrizationOptimalGPU(double *ioU,double *ioV,double error, cpuCompressedMatrixType *initAMatrix,bool directsolver,FILE* logFile);
	//double ParametrizationOptimalGPU_2timeSolves(double *ioU,double *ioV,double error, cpuCompressedMatrixType *initAMatrix,bool directsolver,FILE* logFile);
	double ParametrizationOptimalCPU(double *ioU,double *ioV,double error,int non_zero_element,double *init_sa,unsigned long *init_ija,FILE* logFile);
	#ifdef INTEL_MKL_VERSION
	double ParametrizationOptimalCPU_MKL(double *ioU,double *ioV,double error,int non_zero_element,double *init_sa,unsigned long *init_ija,FILE* logFile);
	#endif
	void SetSigma(double *ipU,double *ipV,double *opSigma ,double gamma = 1.0);
	double GetStretchError(double *ipU,double *ipV);
	double GetStretchError(double *ipU,double *ipV, bool include_boundary, double *opFaceStretch);
	void setSigma(double gamma);
	//IDCutHedge *m_pCutHedgeH;
	//IDCutHedge *m_pCutHedgeT;

	MirrorFace *mirrorFace;
	double *sigmaU;
	double *sigmaV;

	double gammaU;
	double gammaV;

	//void SquareParametrizationGPU(int bottomleftIDvertex,IDList *borderH,IDList *borderT,double total_length_edges, cpuCompressedMatrixType *initA, double *& resultU,double *& resultV,double &resultStretchErr,bool directsolver = false,bool reportlog = true);
	void SquareParametrizationCPU(int bottomleftIDvertex,IDList *borderH,IDList *borderT,double total_length_edges,int non_zero_element, double *init_sa,unsigned long *init_ija , double *& resultU,double *& resultV,double &resultStretchErr,bool reportlog = true);

	//void CircleParametrizationGPU(IDList *borderH,IDList *borderT,double total_length_edges,cpuCompressedMatrixType *initA,  double *& resultU,double *& resultV,double &resultStretchErr,bool directsolver = false,bool reportlog = true);
	void CircleParametrizationCPU(IDList *borderH,IDList *borderT,double total_length_edges,int non_zero_element, double *init_sa,unsigned long *init_ija , double *& resultU,double *& resultV,double &resultStretchErr,bool reportlog = true);

	void StretchAtBoundary(PolarVertex *pIPV, int num_PV,std::vector<double> &op_stretch);
	
};
