#include <stdio.h>
#define Pi 3.141592653589793
void interp3(double **uf, double **uc, int nf,int a) 
{ 
	int ic,iif,jc,jf,nc,i,j; 
	nc=nf/2+1; 
	for (jc=1,jf=1;jc<=nc;jc++,jf+=2)
		for (ic=1;ic<=nc;ic++) 
			uf[2*ic-1][jf]=uc[ic][jc]; 

	for (i=2;i<=nf-1;i++)
		uf[i][1]=0.0;
	for (i=2;i<=nf-1;i++)
		uf[i][nf]=0.0;
	for (j=1;j<=nf;j++)
		uf[1][j]=0.0;
	for (j=1;j<=nf;j++)
		uf[nf][j]=0.0; 
		
	for (jf=3;jf<=nf-2;jf+=2) 
		for (iif=2;iif<=nf-1;iif+=2) 
			uf[iif][jf]=0.5*(uf[iif+1][jf]+uf[iif-1][jf]);


	
	
			
			for (jf=2;jf<=nf-1;jf+=2)
		    for (iif=2;iif<= nf-1;iif++) 
            uf[iif][jf]=0.5*(uf[iif][jf+1]+uf[iif][jf-1]);


	
}