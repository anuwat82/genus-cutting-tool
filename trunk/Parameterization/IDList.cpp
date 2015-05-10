#include "IDList.h"

#ifdef WIN32_MEM_HANDLE
void* IDList::operator new(size_t size)
{
	printf("IDLIST ALLOC NEW\n");
    return  malloc(size);
}

void IDList::operator delete(void* addr)
{
	printf("IDLIST DEALLOC DELETE\n");
    free(addr);
}

#endif //WIN32_MEM_HANDLE
