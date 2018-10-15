#include <stdio.h>
#define Pi 3.141592653589793
#define Up 1.570796326794897
#define Down1 1.570796326794897
#define Down2 1.570796326794897
#define Right 1.570796326794897
#define Left 1.570796326794897
void interp2(double **uf, double **uc, int nf, double m, int a) 
{ 
	int ic,iif,jc,jf,nc,i,j; 
	nc=nf/2+1; 
	for (jc=1,jf=1;jc<=nc;jc++,jf+=2)
		for (ic=1;ic<=nc;ic++) 
			uf[2*ic-1][jf]=uc[ic][jc]; 

	for (i=1;i<=a;i++)
		uf[i][1]=Down1-0.5*m*Pi; 
	
	for (i=1+a;i<=nf;i++)
		uf[i][1]=Down2;
	for (i=1;i<=nf;i++)
		uf[i][nf]=Up-0.5*m*Pi;
	for (j=1;j<=nf;j++)
		uf[nf][j]=Right; 
	for (j=1;j<=nf;j++)
		uf[1][j]=Left-0.5*m*Pi; 


		for (jf=3;jf<=nf-2;jf+=2) 
		for (iif=2;iif<=nf-1;iif+=2) 
			uf[iif][jf]=0.5*(uf[iif+1][jf]+uf[iif-1][jf]);


	
	
			
			for (jf=2;jf<=nf-1;jf+=2)
		    for (iif=2;iif<= nf-1;iif++) 
            uf[iif][jf]=0.5*(uf[iif][jf+1]+uf[iif][jf-1]);



		

}