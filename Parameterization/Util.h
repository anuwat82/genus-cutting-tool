#ifndef __UTIL_HEADER__
#define __UTIL_HEADER__

#include <math.h>

static const double epsilon = 1E-08;

inline void endian_swap16(void* _x)
{
	short* x = (short *)_x;
    *x = ((*x>>8) & 0x00FF) | 
         ((*x<<8) & 0xFF00);
}
inline void endian_swap32(void* _x)
{
	int* x = (int *)_x;
	int y = *x;
    y = ((y>>24)  & 0x000000FF) | 
        ((y<<8)   & 0x00FF0000) |
        ((y>>8)   & 0x0000FF00) |
        ((y<<24)  & 0xFF000000);
	*x = y;
}
// __int64 for MSVC, "long long" for gcc
inline void endian_swap64(void* _x)
{
	__int64* x = (__int64 *)_x;
    *x =((*x>>56) & 0x00000000000000FF) | 
        ((*x<<40) & 0x00FF000000000000) |
        ((*x<<24) & 0x0000FF0000000000) |
        ((*x<<8)  & 0x000000FF00000000) |
        ((*x>>8)  & 0x00000000FF000000) |
        ((*x>>24) & 0x0000000000FF0000) |
        ((*x>>40) & 0x000000000000FF00) |
        ((*x<<56) & 0xFF00000000000000);
}


inline double vector_length(double vec[3])
{
	return sqrt(
				pow(vec[0],2) +
				pow(vec[1],2) +
				pow(vec[2],2)
     			);
}

inline double vector_distance(double *vec1,double *vec2)
{
	return sqrt(
				pow(vec1[0]-vec2[0],2) +
				pow(vec1[1]-vec2[1],2) +
				pow(vec1[2]-vec2[2],2)
     			);
}


inline void vector_normalize(double *src_vec3,double *dst_vec3)
{
	double length = vector_length(src_vec3);

	if (length > 0.0f)
	{
		dst_vec3[0] = src_vec3[0]/length;
		dst_vec3[1] = src_vec3[1]/length;
		dst_vec3[2] = src_vec3[2]/length;
		
	}
	else
	{
		dst_vec3[0]=dst_vec3[1]=dst_vec3[2]=0.0;
	}
}

inline void vector_add(double *vec3_A,double *vec3_B, double *vec3_result)
{
	vec3_result[0] = vec3_A[0] + vec3_B[0];
	vec3_result[1] = vec3_A[1] + vec3_B[1];
	vec3_result[2] = vec3_A[2] + vec3_B[2];
}

inline void vector_scalarMultiply(double *vec3_src,double scalar,double *vec3_dst)
{
	vec3_dst[0] = vec3_src[0]*scalar;
	vec3_dst[1] = vec3_src[1]*scalar;
	vec3_dst[2] = vec3_src[2]*scalar;
}


inline void vector_copy(double *vec3_src,double *vec3_dst)
{
	vec3_dst[0] = vec3_src[0];
	vec3_dst[1] = vec3_src[1];
	vec3_dst[2] = vec3_src[2];
}

inline double vector_dot(double *vec3_A,double *vec3_B)
{
	return ( vec3_A[0]*vec3_B[0] + vec3_A[1]*vec3_B[1] + vec3_A[2]*vec3_B[2] );
}

inline void vector_cross(double *vec3_A,double *vec3_B,double *vec3_result)
{
	vec3_result[0] = ((vec3_A[1])*(vec3_B[2]) - (vec3_B[1])*(vec3_A[2]));
	vec3_result[1] = ((vec3_A[2])*(vec3_B[0]) - (vec3_B[2])*(vec3_A[0]));
	vec3_result[2] = ((vec3_A[0])*(vec3_B[1]) - (vec3_B[0])*(vec3_A[1]));
}

inline bool isPointInsideTriangle(	
									double Px,double Py,
									double x1,double y1,
									double x2,double y2,
									double x3,double y3,
									double *lamda1,
									double *lamda2,
									double *lamda3
								  )
{
	double detT = ((x1-x3)*(y2-y3)) - ((x2-x3)*(y1-y3));
	if (detT == 0)
		return false;
	double _lamda1 = (((y2-y3)*(Px-x3)) + ((x3-x2)*(Py-y3)))/detT;
	double _lamda2 = (((y3-y1)*(Px-x3)) + ((x1-x3)*(Py-y3)))/detT;
	double _lamda3 = 1-_lamda1-_lamda2;

	//double _lamda1 = (((y2-y3)*(Px-x3)) + ((x3-x2)*(Py-y3)));
	//double _lamda2 = (((y3-y1)*(Px-x3)) + ((x1-x3)*(Py-y3)));
	//double _lamda3 = detT-_lamda1-_lamda2;
	
	if (_lamda1 >= -epsilon && _lamda1 <= 1+epsilon &&
		_lamda2 >= -epsilon && _lamda2 <= 1+epsilon &&
		_lamda3 >= -epsilon && _lamda3 <= 1+epsilon)
	/*
	if (_lamda1 >= 0 && _lamda1 <= detT &&
		_lamda2 >= 0 && _lamda2 <= detT &&
		_lamda3 >= 0 && _lamda3 <= detT )
	*/
	{
		
		if (_lamda1 < 0)
			_lamda1 = 0;
		else if (_lamda1 > 1)
			_lamda1 = 1;

		if (_lamda2 < 0)
			_lamda2 = 0;
		else if (_lamda2 > 1)
			_lamda2 = 1;


		if (_lamda3 < 0)
			_lamda3 = 0;
		else if (_lamda3 > 1)
			_lamda3 = 1;
		
		*lamda1 = _lamda1;
		*lamda2 = _lamda2;
		*lamda3 = _lamda3;
		return true;
	}
	else
		return false;
}




#endif //__UTIL_HEADER__