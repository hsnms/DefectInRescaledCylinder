#include "Header2.h"
#include "Header5.h"
double integrat (double **f, int n1, int n2, int n3, int n4, int n)
{
	int i;
	double *k,integ;
	k=dvector(1,n);
	for (i=n1;i<=n2;i++)
	{
		k[i]=trapzd(f[i],n3,n4,n);
	}
	integ=trapzd(k,n1,n2,n);
	free_dvector(k,1,n);
	return integ;
}