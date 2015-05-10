#ifndef __IDLIST_HEADER__
#define __IDLIST_HEADER__

#include <stdio.h>

//#define WIN32_MEM_HANDLE
#ifdef WIN32_MEM_HANDLE

	#include <stdlib.h>
	#include <malloc.h>
#endif 


/*
A fast and simple stretch-minimizing mesh parameterization C++ code
Copyright:(c) Shin Yoshizawa, 2004
E-mail: shin.yoshizawa@mpi-sb.mpg.de
URL: http://www.mpi-sb.mpg.de/~shin
Affiliation: Max-Planck-Institut fuer Informatik: Computer Graphics Group 
 Stuhlsatzenhausweg 85, 66123 Saarbruecken, Germany
 Phone +49 681 9325-408 Fax +49 681 9325-499 

 All right is reserved by Shin Yoshizawa.
This C++ sources are allowed for only primary user of 
research and educational purposes. Don't use secondary: copy, distribution, 
diversion, business purpose, and etc.. 
 */

#define nextN(p) p->next
#define backN(p) p->back

class IDList{
 public:
  int ID;
  IDList *next;
  IDList *back;
  double length;
  IDList(){next=NULL;back=NULL;length=-1.0;}
  IDList(int dv){ID=dv;next=NULL;back=NULL;length=-1.0;}
  virtual ~IDList(){/*ID = -1; next=NULL;back=NULL;length=-1.0;*/}

  #ifdef WIN32_MEM_HANDLE
  public:
        void* operator new(size_t);		
        void operator delete(void*);
  #endif 

 private:
  IDList(const IDList& rhs);
  const IDList &operator=(const IDList& rhs);
   
};





inline IDList* next(IDList *p)
{
	return p->next;
}
inline IDList* back(IDList *p)
{
	return p->back;
}




#endif //__IDLIST_HEADER__