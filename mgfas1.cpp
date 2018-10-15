#include <stdio.h>
#define NPRE /*80*/ 100
#define NPOST /*250*/300
#define ALPHA 0.33
#define NGMAX 15
#include "Header1.h"
#include "Header2.h"
#include "Header3.h"
#include "Header4.h"
void mgfas1(double **u, int n, int maxcyc, double m, int a, double *q)/**/
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
	nn=n/2+1; 
	ngrid=ng-1; 
	irho[ngrid]=dmatrix2(1,nn,1,nn); 
	rstrct2(irho[ngrid],u,nn);  
	while (nn > 3) 
	{ 
		nn=nn/2+1; 
		irho[--ngrid]=dmatrix2(1,nn,1,nn); 
		rstrct2(irho[ngrid],irho[ngrid+1],nn); 
	} 
	nn=3; 
	iu[1]=dmatrix2(1,nn,1,nn); 
	iq[1]=dvector(1,nn);
	for (kk=1,kk2=1;kk<=nn,kk2<=n;kk++,kk2=kk2+(n-1)/(nn-1))
		iq[1][kk]=q[kk2];
	irhs[1]=dmatrix2(1,nn,1,nn); 
	itau[1]=dmatrix2(1,nn,1,nn); 
	itemp[1]=dmatrix2(1,nn,1,nn); 
	slvsm2(n,iu[1],irho[1],m,a,iq[1]);  /*problem*//*printf("NO_1 u[2][2]=%f\t",iu[1][2][2]);*/
	/*copy2(iu[1],irho[1],nn);*/
	free_dmatrix2(irho[1],1,nn,1,nn); 
	ngrid=ng; /*printf("ngrid=%d\n",ngrid);*/
	for (j=2;j<=ngrid;j++) 
	{  /*printf("j=%d\n",j);*/
		nn=2*nn-1; 
		iu[j]=dmatrix2(1,nn,1,nn); 
		iq[j]=dvector(1,nn);
		
		for (kk=1,kk2=1;kk<=nn,kk2<=n;kk++,kk2=kk2+(n-1)/(nn-1))
		iq[j][kk]=q[kk2];
		
		irhs[j]=dmatrix2(1,nn,1,nn); 
		itau[j]=dmatrix2(1,nn,1,nn); 
		itemp[j]=dmatrix2(1,nn,1,nn); 
		interp2(iu[j],iu[j-1],nn,m,a*nn/n); 
		copy2(irhs[j],(j != ngrid ? irho[j] : u),nn); /*printf("irhs[%d][%d][%d]=%f\n",j,nn,nn/2,irhs[j][nn/2][nn/2]);*/
		/*for (k1=1;k1<=nn;k1++)
			for (k2=1;k2<=nn;k2++)
			irhs[j][k1][k2]=0.0;*/
		for (jcycle=1;jcycle<=maxcyc;jcycle++) 
		{ 
			/*printf("jcycle=%d\n",jcycle);*/
			nf=nn; 
			/*printf("j=%d\n",j);*/
			for (jj=j;jj>=2;jj--) 
			{  /*printf("Before jj=%d\n",jj);*/
				for (jpre=1;jpre<=NPRE;jpre++) 
					relax2(iu[jj],irhs[jj],nf,m,a*nf/n,iq[jj]);
				
				lop(itemp[jj],iu[jj],nf,m,a*nf/n,iq[jj]); /**/
			
				

				nf=nf/2+1;
				jm1=jj-1; 
				rstrct2(itemp[jm1],itemp[jj],nf); 
				rstrct2(iu[jm1],iu[jj],nf); 
				lop(itau[jm1],iu[jm1],nf,m,a*nf/n,iq[jm1]); /**/
				matsub(itau[jm1],itemp[jm1],itau[jm1],nf);  
				
				if (jj == j) 
				{trerr=ALPHA*anorm2(itau[jm1],nf); /*printf("trerr=%f\n",trerr);*/}
				rstrct2(irhs[jm1],irhs[jj],nf);  
				matadd(irhs[jm1],itau[jm1],irhs[jm1],nf); 
			}/*printf("irhs[1][2][2]=%f\n",irhs[1][2][2]);*/
			slvsm2(n,iu[1],irhs[1],m,a,iq[1]); /**/
			nf=3; 
			for (jj=2;jj<=j;jj++) 
			{ /*printf("After jj=%d\n",jj);*/
				jm1=jj-1; 
				rstrct2(itemp[jm1],iu[jj],nf); 
				matsub(iu[jm1],itemp[jm1],itemp[jm1],nf);  
				nf=2*nf-1; 
				interp3(itau[jj],itemp[jm1],nf,a*nf/n); 
				matadd(iu[jj],itau[jj],iu[jj],nf); 
			/*for (k1=1;k1<=nf;k1++)
			for (k2=1;k2<=nf;k2++)
				irhs[jj][k1][k2]=0.0;*/
				for (jpost=1;jpost<=NPOST;jpost++)
				{relax2(iu[jj],irhs[jj],nf,m,a*nf/n,iq[jj]);}
				/*printf("j=%d\n",j);*/
			
			} 
			lop(itemp[j],iu[j],nf,m,a*nf/n,iq[j]); /**/
			matsub(itemp[j],irhs[j],itemp[j],nf); 
			res=anorm2(itemp[j],nf); /*printf("res=%f\n",res);*/
			if (res < trerr) break; 
		}
	}
	copy2(u,iu[ngrid],n); 
	for (nn=n,j=ng;j>=1;j--,nn=nn/2+1) 
	{ 
		free_dmatrix2(itemp[j],1,nn,1,nn); 
		free_dmatrix2(itau[j],1,nn,1,nn); 
		free_dmatrix2(irhs[j],1,nn,1,nn); 
		free_dmatrix2(iu[j],1,nn,1,nn); 
		free_dvector(iq[j],1,nn);
		if (j != ng && j != 1) 
			free_dmatrix2(irho[j],1,nn,1,nn); 
	}
}