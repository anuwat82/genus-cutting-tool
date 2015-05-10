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

#include<stdio.h>
#include"IDList.h"
#include"PolarList.h"
#include"IDSet.h"

void IDSet::AppendVF(int myID,IDList *dVTail)
{
    IDList *now = new IDList(myID);
      
    IDList *dummy = dVTail->back;
    now->next = dVTail;
    dVTail->back = now;
    now->back = dummy;
    dummy->next =now;
}


void IDSet::DeleteVF(int pair0,int pair1,IDList *dIHead,IDList *dITail)
{
	IDList *now = dIHead;
	
	while (next(now) != dITail)
	{
		now = now->next;
		if ((now->ID == pair0 && now->next->ID == pair1) ||
			(now->ID == pair1 && now->next->ID == pair0))
		{
			//found!
			now->next->next->back = now->back;
			now->back->next = now->next->next;
			delete now->next;
			delete now;			
			return;
		}
		now = now->next;
	}
      
    
}
int IDSet::SearchI(int dID,IDList *dIHead,IDList *dITail){
    
    IDList *now = dIHead;
    while(now->next!=dITail){
      now = now->next;
      if((now->ID==dID))
	  {
		return 1;
      }
    }
    return 0;
}

IDList * IDSet::GetI(int dID,IDList *dIHead,IDList *dITail)
{
	IDList *now = dIHead;
    while(now->next!=dITail){
      now = now->next;
      if((now->ID==dID))
	  {
		return now;
      }
    }
    return NULL;
}

 
  
 void IDSet::AppendI(int dID,IDList *dIHead,IDList *dITail,int nowID,int *dnum)
{
    if(dID!=nowID)
	{
		if(SearchI(dID,dIHead,dITail)==0)
		{
		  IDList *now = new IDList(dID);
		  IDList *dummy = dITail->back;
		  now->next = dITail;
		  dITail->back = now;
		  now->back = dummy;
		  dummy->next =now;
		  dnum[nowID]++;
		}
		  }else{
	
		  }
}

void IDSet::DeleteI(int dID,IDList *dIHead,IDList *dITail,IDList *dVHead,IDList *dVTail,int nowID,int *dnum )
{
	if(dID!=nowID)
	{
		if (SearchI(dID,dVHead,dVTail)==0)
		{
			IDList *pGetNode = GetI(dID,dIHead,dITail);
			if(pGetNode!=NULL)
			{
				pGetNode->next->back = pGetNode->back;
				pGetNode->back->next = pGetNode->next;
				delete pGetNode;
				if (dnum!=NULL)
					dnum[nowID]--;
				return;
			}
		}
	}
}


void IDSet::AppendPolarI(int dID,PolarList *dITail,double dx,double dy){
  PolarList *now = new PolarList(dID,dx,dy);
  PolarList *dummy = dITail->back;
  now->next = dITail;
  dITail->back = now;
  now->back = dummy;
  dummy->next =now;
}
void IDSet::AppendPolarI(int dID,PolarList *dITail,double dx,double dy,double dz,double cotw){
  PolarList *now = new PolarList(dID,dx,dy,dz,cotw);
  PolarList *dummy = dITail->back;
  now->next = dITail;
  dITail->back = now;
  now->back = dummy;
  dummy->next =now;
}

/*
void IDSet::AppendPolarI(int dID,PolarList *dITail,double dx,double dy,double dz,int nextID,int backID){
  PolarList *now = new PolarList(dID,dx,dy,dz,nextID,backID);
  PolarList *dummy = dITail->back;
  now->next = dITail;
  dITail->back = now;
  now->back = dummy;
  dummy->next =now;
  }*/
/*
void IDSet::AppendPolarIRL(int dID,PolarList *dITail,double dr,double dthta,double dlambda,int tR,int tL,int nextP){
  PolarList *now = new PolarList(dID,dr,dthta,dlambda,tR,tL,nextP);
  PolarList *dummy = dITail->back;
  now->next = dITail;
  dITail->back = now;
  now->back = dummy;
  dummy->next =now;
  


}
*/

void IDSet::DeleteF(int dID,IDList *dIHead,IDList *dITail)
{
	IDList *now = dIHead;
	while(now->next!=dITail)
	{
		now = now->next;
		if((now->ID==dID))
		{
			now->next->back = now->back;
			now->back->next = now->next;
			delete now;
			return;
		}
	}    
}

void IDSet::AppendVFSort(int dID,IDList *dIHead,IDList *dITail)
{
  IDList *now = new IDList(dID);
  IDList *dummy = dITail->back;
  if(dIHead->next!=dITail)
    while(dID<dummy->ID){
      dummy = dummy->back;
      if(dummy==dIHead)break;
    }
  dummy = dummy->next;
  IDList *dummyback=dummy->back;
  
  now->next = dummy;
  dummy->back = now;
  now->back = dummyback;
  dummyback->next =now;
}

void IDSet::AppendISort(int dID,IDList *dIHead,IDList *dITail){
  
  
    
    if(SearchI(dID,dIHead,dITail)==0){
      
      IDList *now = new IDList(dID);
      IDList *dummy = dITail->back;
      if(dIHead->next!=dITail)
	while(dID<dummy->ID){
	  dummy = dummy->back;
	  if(dummy==dIHead)break;
	}
      dummy = dummy->next;
      IDList *dummyback=dummy->back;
      
      now->next = dummy;
      dummy->back = now;
      now->back = dummyback;
      dummyback->next =now;
      
    }
}
/*
void AppendISortLength(int dID,IDList *dIHead,IDList *dITail,double length)
{

}
*/

void IDSet::AppendISort(int dID,IDList *dIHead,IDList *dITail,int nowID,int *dnum)
{
  
	if(dID!=nowID)
	{    
		if(SearchI(dID,dIHead,dITail)==0)
		{
      
			IDList *now = new IDList(dID);
			IDList *dummy = dITail->back;
			if(dIHead->next!=dITail)
			{
				while(dID<dummy->ID)
				{
					dummy = dummy->back;
					if(dummy==dIHead)
						break;
				}
			}
			dummy = dummy->next;
			IDList *dummyback=dummy->back;
      
			now->next = dummy;
			dummy->back = now;
			now->back = dummyback;
			dummyback->next =now;
			dnum[nowID]++;
		}
    
	}
	else
	{
    
	}  
}
void IDSet::AppendIF(int dID,IDList *dIHead,IDList *dITail,int nowID,int *dnum)
{
      
	if(SearchI(dID,dIHead,dITail)==0)
	{
	  IDList *now = new IDList(dID);
	  IDList *dummy = dITail->back;
	  now->next = dITail;
	  dITail->back = now;
	  now->back = dummy;
	  dummy->next =now;
	  dnum[nowID]++;
	}
}


void IDSet::CleanNeighborPolar(PolarList* dHead,PolarList* dTail)
{

	PolarList *dummy = NULL; PolarList *now = NULL;
	if (dHead != NULL && dTail != NULL)
	{		
		now = dHead;
		dummy = now->next;
		while(dummy != dTail)
		{
			now = dummy;
			dummy = now->next;
			delete now;
		}
	}	
	delete dHead;
	delete dTail;


}

void IDSet::CleanNeighbor(IDList* dHead,IDList* dTail)
{
  IDList *dummy = NULL; IDList *now =  NULL;
  if (dHead == dTail)
  {
	  delete dHead;
  }
  else
  {
		now = dHead;
		dummy = now->next;
		while(dummy != dTail)
		{
			now = dummy;
			dummy = now->next;
			delete now;
		}
		delete dHead;
		delete dTail;
  }
}
void IDSet::Clean(IDList **dFHead,IDList **dFTail,int numberSV,int *dneighborN){
  IDList *now=NULL;
  IDList *dummy=NULL;
  int i;
  for(i=0;i<numberSV;i++)
  {
	if(dneighborN[i]!=0)
	{
		now = dFHead[i];
		dummy = now->next;
		while(dummy != dFTail[i])
		{
			now = dummy;
			dummy = now->next;
			delete now;
		}
    }
    
    dFHead[i]->next = dFTail[i];
    dFTail[i]->back = dFHead[i];
    dneighborN[i]=0;
  }
  
}
  
void IDSet::CleanNeighborLL(IDList **dFHead,IDList **dFTail,int numberSV,int *dneighborN){
  IDList *now=NULL;
  IDList *dummy=NULL;
  int i;
  for(i=0;i<numberSV;i++)
  {
    if(dneighborN[i]!=0)
	{
		now = dFHead[i];
		dummy = now->next;
		while(dummy != dFTail[i])
		{
			now = dummy;
			dummy = now->next;
			delete now;
		}
    }
    delete dFHead[i];
    delete dFTail[i]; 
    
  }
  delete [] dneighborN;
  delete [] dFHead;
  delete [] dFTail;
 }
void IDSet::CleanNeighborL(IDList **dFHead,IDList **dFTail,int numberSV){
  IDList *now=NULL;
  IDList *dummy=NULL;
  int i;
  for(i=0;i<numberSV;i++)
  {
	if (dFHead[i])
	{
		now = dFHead[i];
		dummy = now->next;
		while(dummy != dFTail[i])
		{
			now = dummy;
			dummy = now->next;
			delete now;
		}
		delete dFHead[i];
		delete dFTail[i]; 
	}
  }
  delete [] dFHead;
  delete [] dFTail;
  }


void IDSet::CleanNeighborLPolar(PolarList **dFHead,PolarList **dFTail,int numberSV){
  PolarList *now=NULL;
  PolarList *dummy=NULL;
  int i;
  for(i=0;i<numberSV;i++)
  {
    if (dFHead[i])
	{
		now = dFHead[i];
		dummy = now->next;
		while(dummy != dFTail[i])
		{
			now = dummy;
			dummy = now->next;
			delete now;
		}
		delete dFHead[i];
		delete dFTail[i]; 
	}    
  }
  delete [] dFHead;
  delete [] dFTail;
 }





