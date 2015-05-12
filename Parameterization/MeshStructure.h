#pragma once
#include <stdio.h>





typedef struct Vertex 
{

	double x;
	double y;
	double z;
	Vertex()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}

	Vertex(double _x,double _y, double _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}

	void Set(double *xyz)	
	{
		x = xyz[0];
		y = xyz[1];
		z = xyz[2];
	}
	
} Vertex;


typedef struct PolarVertex
{
	Vertex *p_vertex;
	double u;
	double v;
	PolarVertex()
	{
		p_vertex = NULL;
		u = 0.0;
		v = 0.0;
	}
} PolarVertex;

typedef struct Face 
{
	unsigned char nverts;    /* number of vertex indices in list */	
	int *faceInfo_ptr;       /* pointer to the first array of this face information*/

	double normalVector[3];
	Face()
	{	
		nverts = 0;
		faceInfo_ptr = NULL;
		normalVector[0]=normalVector[1]=normalVector[2] = 0.0f;
	}
} Face;