#include "GPUSolver.cuh"
#include "MyParameterization.h"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

//#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/matrix_proxy.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/linalg/direct_solve.hpp"
#include "viennacl/linalg/lu.hpp"
#include "viennacl/linalg/ichol.hpp"

//Parallel Patterns Library
#include <ppl.h>
#include <concrt.h>


typedef double ScalarType;
//typedef boost::numeric::ublas::compressed_matrix<ScalarType>        cpuCompressedMatrixType;
typedef boost::numeric::ublas::compressed_matrix<ScalarType>        cpuCompressedMatrixType;
typedef boost::numeric::ublas::matrix<ScalarType>					cpuDenseMatrixType;
typedef boost::numeric::ublas::vector<ScalarType>					cpuVectorType;
typedef viennacl::compressed_matrix<ScalarType>						gpuCompressedMatrixType;
typedef viennacl::matrix<ScalarType>								gpuDenseMatrixType;
typedef viennacl::vector<ScalarType>								gpuVectorType;

using namespace boost::numeric;
void MyParameterization::ParametrizationOptimalGPU(double error,FILE* logFile)
{

	//int numberV = 100;
	cpuCompressedMatrixType cpuAmatrix(numberV*2,numberV*2);
	cpuVectorType cpuBvector(numberV*2);
	cpuVectorType cpuXvector(numberV*2);
	
	gpuCompressedMatrixType gpuAmatrix(numberV*2,numberV*2);
	gpuVectorType gpuBvector(numberV*2);
	gpuVectorType gpuXvector(numberV*2);
	
	setFloaterC();
	SortIndexP();

	PolarList *nowp = NULL;
	for (int i = 0; i < numberV ; i++)
	{
		if(boundary[i]!=1)
		{
			cpuBvector(i) = 0.0f;
			cpuBvector(i+numberV) = 0.0f;

			nowp = PHead[i];
			
			cpuAmatrix(i,i) = 1.0f;
			cpuAmatrix(numberV +i,numberV +i) = 1.0f;
			while(nextN(nowp)!=PTail[i])
			{
				nowp = nextN(nowp);						
				cpuAmatrix(i,nowp->ID) = -nowp->lambda;		
				cpuAmatrix(numberV +i,numberV + nowp->ID) = -nowp->lambda;
			}			
		}
		else
		{			
			//constraint 
			cpuBvector(i) = pU[i];
			cpuBvector(i+numberV) = pV[i];

			cpuAmatrix(i,i) = 1.0f;
			cpuAmatrix(numberV + i,numberV +i) = 1.0f;
		}
	}
			
	
	viennacl::copy(cpuAmatrix, gpuAmatrix);
	viennacl::copy(cpuBvector, gpuBvector);
	
			
	viennacl::linalg::ilu0_tag ilu0_config;
	viennacl::linalg::ilu0_precond< gpuCompressedMatrixType > vcl_ilut(gpuAmatrix,ilu0_config);
	clock_t calTime = 0;
	
	gpuXvector = viennacl::linalg::solve(gpuAmatrix,gpuBvector,viennacl::linalg::bicgstab_tag(error),vcl_ilut);
	calTime = clock() - calTime; 
	viennacl::copy(gpuXvector, cpuXvector);

	for(int i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			pU[i] = cpuXvector(i);
			pV[i] = cpuXvector(numberV + i);
					
		}
	}
	
	
	
	
	double previous_l2= 0;
	double candidate_l2 = this->getCurrentE();
	this->resultStretch = candidate_l2;
	double *sigsum = new double[numberV];
	double *prevU = new double[numberV];
	double *prevV = new double[numberV];
	do
	{
		previous_l2 = candidate_l2;
		memcpy(prevU,pU,sizeof(double)*numberV);
		memcpy(prevV,pV,sizeof(double)*numberV);
		setSigmaZero();
		for(int i=0;i<numberV;i++)
		{      
			if(boundary[i]!=1)
			{
				sigsum[i]=0.0;
				nowp = PHead[i];
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
					nowp->old_lambda = nowp->lambda;
					nowp->lambda /= sigma[nowp->ID];
					sigsum[i] += nowp->lambda;
				}

				cpuAmatrix(i,i) = 1.0f;
				cpuAmatrix(i+numberV,i+numberV) = 1.0f;
				nowp = PHead[i];
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
					cpuAmatrix(i,nowp->ID) = -nowp->lambda/sigsum[i];	
					cpuAmatrix(i+numberV,nowp->ID+numberV) = -nowp->lambda/sigsum[i];	
				}
			}
			else
			{
				cpuAmatrix(i,i) = 1.0f;
				cpuAmatrix(i+numberV,i+numberV) = 1.0f;
			}
		}

		viennacl::copy(cpuAmatrix, gpuAmatrix);
		viennacl::copy(cpuBvector, gpuBvector);

		viennacl::linalg::ilu0_tag ilu0_config;
		viennacl::linalg::ilu0_precond< gpuCompressedMatrixType > vcl_ilut(gpuAmatrix,ilu0_config);
		clock_t calTime = 0;
		calTime = clock();
		gpuXvector = viennacl::linalg::solve(gpuAmatrix,gpuBvector,viennacl::linalg::bicgstab_tag(error),vcl_ilut);
		calTime = clock() - calTime; 
		viennacl::copy(gpuXvector, cpuXvector);

		for(int i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				
				pU[i] = cpuXvector(i);
				pV[i] = cpuXvector(numberV + i);
					
			}
		}		
		candidate_l2 = this->getCurrentE();
	}
	while (previous_l2 > candidate_l2);

	memcpy(pU,prevU,sizeof(double)*numberV);
	memcpy(pV,prevV,sizeof(double)*numberV);

	delete [] prevU;
	delete [] prevV;
	this->resultStretch = previous_l2;
}


double MyParameterization::ParametrizationOptimalGPU(double *ioU,double *ioV,double error, cpuCompressedMatrixType *initAMatrix,bool directsolver,FILE* logFile)
{
	cpuCompressedMatrixType cpuAmatrix(numberV*2,numberV*2);
	if (initAMatrix != NULL && 
		initAMatrix->size1() == initAMatrix->size2() &&
		initAMatrix->size1() == numberV)
	{
		ublas::project(cpuAmatrix,ublas::range(0, numberV),ublas::range(0, numberV)) = *initAMatrix;
		ublas::project(cpuAmatrix,ublas::range(numberV, 2*numberV),ublas::range(numberV, 2*numberV)) = *initAMatrix;
		//ublas::subrange(cpuAmatrix,0,numberV,0,numberV) = *initAMatrix;
		//ublas::subrange(cpuAmatrix,numberV,numberV,numberV,numberV) = *initAMatrix;
	}
	
	cpuVectorType cpuBvector(numberV*2);
	cpuVectorType cpuXvector(numberV*2);
	
	gpuCompressedMatrixType gpuAmatrix(numberV*2,numberV*2);
	gpuVectorType gpuBvector(numberV*2);
	gpuVectorType gpuXvector(numberV*2);

	PolarList *nowp = NULL;
	for (int i = 0; i < numberV ; i++)
	{
		if(boundary[i]!=1)
		{
			cpuBvector(i) = 0.0f;
			cpuBvector(i+numberV) = 0.0f;
			if (initAMatrix == NULL)
			{
				nowp = PHead[i];			
				cpuAmatrix(i,i) = 1.0f;
				cpuAmatrix(numberV +i,numberV +i) = 1.0f;
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);						
					cpuAmatrix(i,nowp->ID) = -nowp->lambda;		
					cpuAmatrix(numberV +i,numberV + nowp->ID) = -nowp->lambda;
				}
			}
		}
		else
		{			
			//constraint 
			cpuBvector(i) = ioU[i];
			cpuBvector(i+numberV) = ioV[i];
			if (initAMatrix == NULL)
			{
				cpuAmatrix(i,i) = 1.0f;
				cpuAmatrix(numberV + i,numberV +i) = 1.0f;
			}
		}
	}
			
	
	
	if (!directsolver)
	{
		viennacl::copy(cpuAmatrix, gpuAmatrix);
		viennacl::copy(cpuBvector, gpuBvector);				
		viennacl::linalg::ilu0_tag ilu0_config;
		viennacl::linalg::ilu0_precond< gpuCompressedMatrixType > vcl_ilut(gpuAmatrix,ilu0_config);
		gpuXvector = viennacl::linalg::solve(gpuAmatrix,gpuBvector,viennacl::linalg::bicgstab_tag(error),vcl_ilut);
		//gpuXvector = viennacl::linalg::solve(gpuAmatrix,gpuBvector,viennacl::linalg::bicgstab_tag(error));
		viennacl::copy(gpuXvector, cpuXvector);
		
	}
	else
	{
		cpuDenseMatrixType denMatA_cpu = cpuAmatrix;
		gpuDenseMatrixType denMatA_gpu ;
		viennacl::copy(denMatA_cpu, denMatA_gpu);	
		viennacl::copy(cpuBvector, gpuBvector);	
	
		//boost::numeric::ublas::permutation_matrix<size_t> pm ( denMatA_cpu.size1() );
		//boost::numeric::ublas::lu_factorize( denMatA_cpu, pm );
		//boost::numeric::ublas::lu_substitute( denMatA_cpu, pm, cpuXvector );

		viennacl::linalg::lu_factorize(denMatA_gpu);
		viennacl::linalg::lu_substitute(denMatA_gpu,gpuBvector); //gpuBvector  got result
		viennacl::copy(gpuBvector, cpuXvector);
	}


	
	

	for(int i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			ioU[i] = cpuXvector(i);
			ioV[i] = cpuXvector(numberV + i);
					
		}
	}
	
	
	double *_sigma  = new double[numberV];	
	double candidate_l2 = GetStretchError(ioU,ioV);
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
		for(int i=0;i<numberV;i++)
		{      
			if(boundary[i]!=1)
			{
				sigsum[i]=0.0;
				nowp = PHead[i];
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
					//nowp->old_lambda = nowp->lambda;
					//nowp->lambda /= sigma[nowp->ID];
					
					sigsum[i] += (-cpuAmatrix(i,nowp->ID)/_sigma[nowp->ID]);
				}

				cpuAmatrix(i,i) = 1.0f;
				cpuAmatrix(i+numberV,i+numberV) = 1.0f;
				nowp = PHead[i];
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
					double newLamda = ((cpuAmatrix(i,nowp->ID)/_sigma[nowp->ID]))/sigsum[i];
					cpuAmatrix(i,nowp->ID) = newLamda;	
					cpuAmatrix(i+numberV,nowp->ID+numberV) = newLamda;	
				}
			}
			else
			{
				cpuAmatrix(i,i) = 1.0f;
				cpuAmatrix(i+numberV,i+numberV) = 1.0f;
			}
		}

		if (!directsolver)
		{
			viennacl::copy(cpuAmatrix, gpuAmatrix);
			viennacl::copy(cpuBvector, gpuBvector);				
			viennacl::linalg::ilu0_tag ilu0_config;
			viennacl::linalg::ilu0_precond< gpuCompressedMatrixType > vcl_ilut(gpuAmatrix,ilu0_config);
			gpuXvector = viennacl::linalg::solve(gpuAmatrix,gpuBvector,viennacl::linalg::bicgstab_tag(error),vcl_ilut);
			//gpuXvector = viennacl::linalg::solve(gpuAmatrix,gpuBvector,viennacl::linalg::bicgstab_tag(error));
			viennacl::copy(gpuXvector, cpuXvector);
		
		}
		else
		{
			cpuDenseMatrixType denMatA_cpu = cpuAmatrix;
			gpuDenseMatrixType denMatA_gpu ;
			viennacl::copy(denMatA_cpu, denMatA_gpu);	
			viennacl::copy(cpuBvector, gpuBvector);	
			viennacl::linalg::lu_factorize(denMatA_gpu);
			viennacl::linalg::lu_substitute(denMatA_gpu,gpuBvector); //gpuBvector  got result
			viennacl::copy(gpuBvector, cpuXvector);
		}

		for(int i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				
				ioU[i] = cpuXvector(i);
				ioV[i] = cpuXvector(numberV + i);
					
			}
		}		
		candidate_l2 = GetStretchError(ioU,ioV);
	}
	while (previous_l2 > candidate_l2);

	memcpy(ioU,prevU,sizeof(double)*numberV);
	memcpy(ioV,prevV,sizeof(double)*numberV);

	delete [] _sigma;
	delete [] sigsum;
	delete [] prevU;
	delete [] prevV;
	return previous_l2;

}

double MyParameterization::ParametrizationOptimalGPU_2timeSolves(double *ioU,double *ioV,double error, cpuCompressedMatrixType *initAMatrix,bool directsolver,FILE* logFile)
{
	bool hasInitMatrix = false;
	cpuCompressedMatrixType cpuAmatrix(numberV,numberV);
	if (initAMatrix != NULL && 
		initAMatrix->size1() == initAMatrix->size2() &&
		initAMatrix->size1() == numberV)
	{
		cpuAmatrix = *initAMatrix;
		hasInitMatrix = true;
	}
	gpuCompressedMatrixType gpuAmatrix(numberV,numberV);	


	cpuVectorType cpuB1vector(numberV);
	cpuVectorType cpuB2vector(numberV);
	cpuVectorType cpuX1vector(numberV);	
	cpuVectorType cpuX2vector(numberV);	

	gpuVectorType gpuB1vector(numberV);
	gpuVectorType gpuB2vector(numberV);
	gpuVectorType gpuX1vector(numberV);	
	gpuVectorType gpuX2vector(numberV);	
	

	PolarList *nowp = NULL;
	for (int i = 0; i < numberV ; i++)
	{
		if(boundary[i]!=1)
		{
			cpuB1vector(i) = 0.0f;
			cpuB2vector(i) = 0.0f;
			if (!hasInitMatrix)
			{
				nowp = PHead[i];			
				cpuAmatrix(i,i) = 1.0f;
				
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);						
					cpuAmatrix(i,nowp->ID) = -nowp->lambda;	
				}
			}
		}
		else
		{			
			//constraint 
			cpuB1vector(i) = ioU[i];
			cpuB2vector(i) = ioV[i];			
			if (!hasInitMatrix)
			{
				cpuAmatrix(i,i) = 1.0f;
			}
		}
	}

	
	
	if (!directsolver)
	{
		viennacl::copy(cpuAmatrix, gpuAmatrix);
		viennacl::copy(cpuB1vector, gpuB1vector);
		viennacl::copy(cpuB2vector, gpuB2vector);
		
		viennacl::linalg::ilu0_tag ilu0_config;
		viennacl::linalg::ilu0_precond< gpuCompressedMatrixType > vcl_ilut(gpuAmatrix,ilu0_config);
		//viennacl::linalg::ichol0_tag ichol0_config;
		//viennacl::linalg::ichol0_precond< gpuCompressedMatrixType > vcl_ichol(gpuAmatrix,ichol0_config);
		
		gpuX1vector = viennacl::linalg::solve(gpuAmatrix,gpuB1vector,viennacl::linalg::bicgstab_tag(error),vcl_ilut);
		gpuX2vector = viennacl::linalg::solve(gpuAmatrix,gpuB2vector,viennacl::linalg::bicgstab_tag(error),vcl_ilut);
		//gpuX1vector = viennacl::linalg::solve(gpuAmatrix,gpuB1vector,viennacl::linalg::gmres_tag(error),vcl_ichol);
		//gpuX2vector = viennacl::linalg::solve(gpuAmatrix,gpuB2vector,viennacl::linalg::gmres_tag(error),vcl_ichol);
		
		viennacl::copy(gpuX1vector, cpuX1vector);
		viennacl::copy(gpuX2vector, cpuX2vector);
		
	}
	else
	{
		
		cpuDenseMatrixType denMatA_cpu( cpuAmatrix);
		gpuDenseMatrixType denMatA_gpu(numberV,numberV);
		viennacl::copy(denMatA_cpu, denMatA_gpu);
		viennacl::copy(cpuB1vector, gpuB1vector);
		viennacl::copy(cpuB2vector, gpuB2vector);

		viennacl::linalg::lu_factorize(denMatA_gpu);
		viennacl::linalg::lu_substitute(denMatA_gpu,gpuB1vector); //gpuB1vector  got result
		viennacl::linalg::lu_substitute(denMatA_gpu,gpuB2vector); //gpuB2vector  got result
		
		viennacl::copy(gpuB1vector, cpuX1vector);
		viennacl::copy(gpuB2vector, cpuX2vector);
		
	}


	
#if 1

	for(int i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			ioU[i] = cpuX1vector(i);
			ioV[i] = cpuX2vector(i);
		}
	}
	
	
	double *_sigma  = new double[numberV];	
	double candidate_l2 = GetStretchError(ioU,ioV);
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
		for(int i=0;i<numberV;i++)
		{      
			if(boundary[i]!=1)
			{
				sigsum[i]=0.0;
				nowp = PHead[i];
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
					//nowp->old_lambda = nowp->lambda;
					//nowp->lambda /= sigma[nowp->ID];
					
					sigsum[i] += (-cpuAmatrix(i,nowp->ID)/_sigma[nowp->ID]);
				}

				cpuAmatrix(i,i) = 1.0f;
				
				nowp = PHead[i];
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
					double newLamda = ((cpuAmatrix(i,nowp->ID)/_sigma[nowp->ID]))/sigsum[i];
					cpuAmatrix(i,nowp->ID) = newLamda;
				}
			}
			else
			{
				cpuAmatrix(i,i) = 1.0f;
			}
		}

		if (!directsolver)
		{
			viennacl::copy(cpuAmatrix, gpuAmatrix);
			viennacl::copy(cpuB1vector, gpuB1vector);
			viennacl::copy(cpuB2vector, gpuB2vector);
		
			viennacl::linalg::ilu0_tag ilu0_config;
			viennacl::linalg::ilu0_precond< gpuCompressedMatrixType > vcl_ilut(gpuAmatrix,ilu0_config);
			gpuX1vector = viennacl::linalg::solve(gpuAmatrix,gpuB1vector,viennacl::linalg::bicgstab_tag(error),vcl_ilut);
			gpuX2vector = viennacl::linalg::solve(gpuAmatrix,gpuB2vector,viennacl::linalg::bicgstab_tag(error),vcl_ilut);
		
			viennacl::copy(gpuX1vector, cpuX1vector);
			viennacl::copy(gpuX2vector, cpuX2vector);
		
		}
		else
		{
			cpuDenseMatrixType denMatA_cpu( cpuAmatrix);
			gpuDenseMatrixType denMatA_gpu(numberV,numberV);
			viennacl::copy(denMatA_cpu, denMatA_gpu);
			viennacl::copy(cpuB1vector, gpuB1vector);
			viennacl::copy(cpuB2vector, gpuB2vector);

			viennacl::linalg::lu_factorize(denMatA_gpu);
			viennacl::linalg::lu_substitute(denMatA_gpu,gpuB1vector); //gpuB1vector  got result
			viennacl::linalg::lu_substitute(denMatA_gpu,gpuB2vector); //gpuB2vector  got result
		
			viennacl::copy(gpuB1vector, cpuX1vector);
			viennacl::copy(gpuB2vector, cpuX2vector);
		
			
		}

		for(int i=0;i<numberV;i++)
		{
			if(boundary[i]!=1)
			{
				ioU[i] = cpuX1vector(i);
				ioV[i] = cpuX2vector(i);
			}
		}		
		candidate_l2 = GetStretchError(ioU,ioV);
	}
	while (previous_l2 > candidate_l2);

	memcpy(ioU,prevU,sizeof(double)*numberV);
	memcpy(ioV,prevV,sizeof(double)*numberV);

	delete [] _sigma;
	delete [] sigsum;
	delete [] prevU;
	delete [] prevV;
	return previous_l2;
	#else
	return 0;
#endif
	
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

void ParametrizationOptimalGPU(Polyhedron *poly,FILE* logFile)
{
	int numberV = poly->numberV;
	//int numberV = 100;
	cpuCompressedMatrixType cpuAmatrix(numberV*2,numberV*2);
	cpuVectorType cpuBvector(numberV*2);
	cpuVectorType cpuXvector(numberV*2);
	
	gpuCompressedMatrixType gpuAmatrix(numberV*2,numberV*2);
	gpuVectorType gpuBvector(numberV*2);
	gpuVectorType gpuXvector(numberV*2);
	
	
	
	
	for(int i=0;i<numberV;i++)
	{		
		if(poly->boundary[i]!=1)
		{			
			cpuBvector(i) = 0.0f;
			cpuBvector(i+numberV) = 0.0f;			
		}
		else
		{
			//constraint 
			cpuBvector(i) = poly->pU[i];
			cpuBvector(i+numberV) = poly->pV[i];
		}
	}
	//setFloaterC();
	//poly->SortIndexP();

	PolarList *nowp = NULL;

	//for U
	nowp = NULL;
	for (int i = 0; i < numberV ; i++)
	{
		if(poly->boundary[i]!=1)
		{
			nowp = poly->PHead[i];
			
			cpuAmatrix(i,i) = 1.0f;
			while(nextN(nowp)!=poly->PTail[i])
			{
				nowp = nextN(nowp);						
				cpuAmatrix(i,nowp->ID) = -nowp->lambda;				
			}			
		}
		else
		{			
			cpuAmatrix(i,i) = 1.0f;
		}
	}
			
	//for V
	nowp = NULL;
	for (int i = 0; i < numberV ; i++)
	{
		if(poly->boundary[i]!=1)
		{
			nowp = poly->PHead[i];			
			cpuAmatrix(numberV +i,numberV +i) = 1.0f;			
			while(nextN(nowp)!=poly->PTail[i])
			{
				nowp = nextN(nowp);
				cpuAmatrix(numberV +i,numberV + nowp->ID) = -nowp->lambda;
				
			}			
		}
		else
		{
			cpuAmatrix(numberV + i,numberV +i) = 1.0f;				
		}
	}
	viennacl::copy(cpuAmatrix, gpuAmatrix);
	viennacl::copy(cpuBvector, gpuBvector);

	gpuXvector = viennacl::linalg::solve(gpuAmatrix,gpuBvector,viennacl::linalg::bicgstab_tag());
	viennacl::copy(gpuXvector, cpuXvector);

	for(int i=0;i<numberV;i++)
	{
		if(poly->boundary[i]!=1)
		{
			poly->pU[i] = cpuXvector(i);
			poly->pV[i] = cpuXvector(numberV + i);
					
		}
	}
	//double initialstrech = poly->getCurrentE();
	
	

#if 0

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

	setFloaterC();
	/*
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
	*/
    
  
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
		UaXY[i+1] = pU[i];
		UaXY[i+numberV+1] = pV[i];
	}
	mybcg->linbcg(((unsigned long)(2*(numberV))),vecb,UaXY,1,error,itenum,&iter,&linerr);
  
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			pU[i] = UaXY[i+1];
			pV[i] = UaXY[i+numberV+1];
		}
	}
  
	// Re-solving linear system
	double initialstrech=0.0;
	double currentstrech=0.0;
	Point2d *prevU = new Point2d [numberV];  

	initialstrech = getCurrentE();
	int kk = 0;  

	printf("U%d  STRETCH: %f\n",kk,initialstrech); //u0
	if (logFile)
	{
		fprintf(logFile,"U%d  STRETCH: %f\n",kk,initialstrech); //u0
	}
	for(kk=0;kk<itenum;kk++)
	{
		setSigmaZero();
		for(i=0;i<numberV;i++)
		{
      
			if(boundary[i]!=1)
			{
				sigsum[i]=0.0;
				nowp = PHead[i];
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
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
	
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
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
				while(nextN(nowp)!=PTail[i])
				{
					nowp = nextN(nowp);
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
		printf("currentstrech U%d= %lf\n",kk+1,currentstrech);
		printf("U%d  STRETCH: %f\n",kk+1,currentstrech); 

		if (logFile)
		{
			fprintf(logFile,"U%d  STRETCH: %f\n",kk+1,currentstrech); 
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

			resultStretch = initialstrech;
			break;		
		}
		else
		{
			initialstrech = currentstrech;
		}
	}
  
	level = kk;
 
 
	printf("STL2 (U%d) error = %lf\n", level ,resultStretch);
	delete mybcg;
	/*
	for(i=0;i<numberV;i++)
	{
		if(boundary[i]!=1)
		{
			IDtool->CleanNeighborPolar(PHead[i],PTail[i]);
      
			PHead[i] = new PolarList();
			PTail[i] = new PolarList();
			PHead[i]->next = PTail[i];
			PTail[i]->back = PHead[i];
		}
	}
	*/

	delete [] prevU;
	delete [] vecb;
	delete [] UaXY;
	delete [] sigsum;
#endif
}




double    MyParameterization::PARAM_PARALLEL_GPU(	PolarVertex *pIPV,
												int num_PV,
												FILE* logFile)
{
	

	
	iteNum = (pow((double)((numberV/20000) + 1),2)) *2000;	
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
	/*
	cpuCompressedMatrixType initialAMatrix(numberV*2,numberV*2);
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
	*/
	
	cpuCompressedMatrixType initialAMatrix(numberV,numberV);
	for (int i = 0; i < numberV ; i++)
	{
		if(boundary[i]!=1)
		{
			
			PolarList *nowp = PHead[i];			
			initialAMatrix(i,i) = 1.0f;			
			while(nextN(nowp)!=PTail[i])
			{
				nowp = nextN(nowp);						
				initialAMatrix(i,nowp->ID) = -nowp->lambda;						
			}			
		}
		else
		{			
			//constraint
			initialAMatrix(i,i) = 1.0f;
		}
	}
	
	double *bestU = NULL;
	double *bestV = NULL;
	std::vector<double *>_resultU(bottomLeftVertexIDList.size(),NULL);
	std::vector<double *>_resultV(bottomLeftVertexIDList.size(),NULL);
	std::vector<double>_resultError(bottomLeftVertexIDList.size(),0.0);
	
	unsigned int n = bottomLeftVertexIDList.size(); 
#if 1
	
	concurrency::SchedulerPolicy oldpolicy = concurrency::CurrentScheduler::GetPolicy();	
	concurrency::SchedulerPolicy policy(oldpolicy);
	if (policy.GetPolicyValue(concurrency::MaxConcurrency) > 10)
		policy.SetConcurrencyLimits(1,10);	
	
	concurrency::CurrentScheduler::Create(policy);	
	concurrency::parallel_for(0u, n, [&_resultU,&_resultV,&_resultError,&initialAMatrix,bottomLeftVertexIDList,BpointH,BpointT,tlength,this](unsigned int i)  
	{
		SquareParametrizationGPU(bottomLeftVertexIDList[i], BpointH,BpointT, tlength,&initialAMatrix,_resultU[i],_resultV[i],_resultError[i],true,false);		
	});
	concurrency::CurrentScheduler::Detach();
	
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
	
#else
	
	for (int i = 0 ; i < bottomLeftVertexIDList.size(); i++)
	{
		double *resultU = NULL;
		double *resultV = NULL;
		double resultErr  = 0;
		SquareParametrizationGPU (bottomLeftVertexIDList[i], BpointH,BpointT, tlength,&initialAMatrix,resultU,resultV,resultErr,false,false);
		if (resultErr < bestStretch )
		{
			if (bestU)
				delete [] bestU;
			if (bestV)
				delete [] bestV;
			bestU = resultU;
			bestV = resultV;
			
			bestStretch = resultErr;
			best_startID = i;
		}
		else
		{
			
			delete [] resultU;
			delete [] resultV;
		}

		if (resultErr > worstStretch)
		{
			worstStretch = resultErr;
			worst_startID = i;
		}		
		
	}
	
#endif	
	for (int i=0;i<numberV;i++)
	{				
		pIPV[i].u = bestU[i];
		pIPV[i].v = bestV[i];
	}
	double *resultCircleU = NULL;
	double *resultCircleV = NULL;
	double resultCircleErr = 0;
	CircleParametrizationGPU(BpointH,BpointT, tlength,&initialAMatrix,resultCircleU,resultCircleV,resultCircleErr,false,false);


	if (bestU)
		delete [] bestU;
	if (bestV)
		delete [] bestV;

	if (resultCircleU)
		delete [] resultCircleU;
	if (resultCircleV)
		delete [] resultCircleV;
	
	printf("=== Find best stretch  of TEST%d  (best corner at %f) ===\n",best_startID,bestStretch);
	if (logFile)
		fprintf(logFile,"=== Find best stretch of TEST%d (best corner ERR = %f)===\n",best_startID,bestStretch);
	
	printf("=== Find worst stretch  of TEST%d  (worst corner at %f) ===\n",worst_startID,worstStretch);
	if (logFile)
		fprintf(logFile,"=== Find worst stretch  of TEST%d  (worst corner ERR = %f) ===\n",worst_startID,worstStretch);

	printf("=== Circular Parameterization  stretch = %f ===\n",resultCircleErr);
	if (logFile)
		fprintf(logFile,"=== Circular Parameterization  stretch = %f ===\n",resultCircleErr);

		
	IDtool->CleanNeighbor(BpointH,BpointT);
	
	this->resultStretch = bestStretch;
	return bestStretch;	
}