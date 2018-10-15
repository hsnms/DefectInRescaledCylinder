#include <stdlib.h>
#include "Header5.h"
#define NR_END 1
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror3("allocation failure in ivector()");
	return v-nl+NR_END;
}