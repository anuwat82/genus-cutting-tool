#include "API.h"
#include <memory>
#include<stdio.h>
#include<math.h>
#include"Point3d.h"
#include"Point2d.h"
#include"IDList.h"
#include"PointTool.h"
#include"PolarList.h"
#include"IDSet.h"
#include"PCBCGSolver.h"
#include"Polyhedron.h"

double NaturalParameterizationOptimal(double *vertice, int number_vertice , int *triangleIndices, int number_triangle, double *opU, double *opV)
{

  original::Polyhedron *mymesh = new original::Polyhedron();

  mymesh->readmesh(vertice, number_vertice , triangleIndices, number_triangle);
  mymesh->boundarytype = 2;  
  mymesh->boundarysigma = 1;
  mymesh->weighttype = 0;
  mymesh->smooth  = 1;
  

  mymesh->param();

  memcpy_s(opU, sizeof(double)* number_vertice, mymesh->pU, sizeof(double)* number_vertice);
  memcpy_s(opV, sizeof(double)* number_vertice, mymesh->pV, sizeof(double)* number_vertice);
  double L2error = mymesh->L2Error;
  delete mymesh;
  return L2error;
}