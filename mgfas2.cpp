#include <stdio.h>
#define NPRE 30 
#define NPOST 1530
#define ALPHA 0.33
#define NGMAX 15
#include "Header1.h"
#include "Header2.h"
#include "Header3.h"
#include "Header4.h"
void mgfas2(double **u, int n, int maxcyc, double m, int a, double *q)
{
	int j,jcycle,jj,jm1,jpost,jpre,nf,ng=0,ngrid,nn,k1,k2,kk,kk2,ff,ff2; 
	double **irho[NGMAX+1],**irhs[NGMAX+1],**itau[NGMAX+1], **itemp[NGMAX+1],**iu[NGMAX+1],*iq[NGMAX+1]; 
	double res,trerr;
	double rtsafe( double H, double h,int n,double x1, double x2, double xcc,double rrhs, double **u, double *q);  
	nn=n; 
	while (nn >>= 1) ng++; 
	if (n != 1+(1L << ng)) 
		nrerror2("n-1 must be a power of 2 in mgfas."); 
	if (ng > NGMAX) 
		nrerror2("increase NGMAX in mglin."); 
	ngrid=ng; 
 j=ngrid;/*printf("j=%d\n",j);*/

	
		
		irhs[j]=dmatrix2(1,n,1,n); 

		
		for (k1=1;k1<=n;k1++)
			for (k2=1;k2<=n;k2++)
				irhs[j][k1][k2]=0.0;

				for (jpost=1;jpost<=NPOST;jpost++)
				{relax2(u,irhs[j],n,m,a,q);/*printf("iu[%d][%d-1][%d-1]=%f\n",jj,nf,nf,iu[jj][nf-1][nf-1]);*/}
				
		
		

	
		
		free_dmatrix2(irhs[j],1,n,1,n); 
	
	
	
}