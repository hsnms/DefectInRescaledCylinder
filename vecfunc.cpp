#include <stdlib.h>
#include <math.h>
#include "Header2.h"
/*#define U 0.0;
#define D 0.0;*/
#define R 1.570796326794897;
/*#define L -1.57;*/
void vecfunc(int n, int k, float *x, float *f, double m, int a)
{
	int i;
	/*f[2]=-(20.0+4.0*16)*x[2]+8.0*x[1]-sin(2.0*x[2])*//*+2.0*10*3.14;*/
	/*f[1]=-(6.0+16.0)*x[1]-sin(2.0*x[1])+6.0*x[2]+0.5*3.14*10;*/
	/*f[1]=2.0*2.0*1.57*m*a/n-2.0*1.57*m-(6.0+4*16*0.25)*x[1]-sin(2.0*x[1]);*/
	for (i=/*2*/1;i<=k;i++)
		f[i]=0.0;
}