#include <stdlib.h>
#include "Header5.h"
#define NR_END 1
float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror3("allocation failure in vector()");
	return v-nl+NR_END;
}