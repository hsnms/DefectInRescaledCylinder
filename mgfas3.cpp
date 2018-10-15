#include <stdio.h>
#include <math.h>
#define NPRE 30 
#define NPOST /*1300*//*800*/1300
#define ALPHA 0.33
#define NGMAX 15
#define Pi 3.141592653589793
#include "Header1.h"
#include "Header2.h"
#include "Header3.h"
#include "Header4.h"
void mgfas3(double **u, int n, int maxcyc, double m, int a,int b, double *q, double *x, double *y, double R0)
{
	int j,i,jpost,lsw,jsw,ipass; 
	double h,res,v1,v2,v3,v4,v5,v6,v7,v8,v9,H,h2i,foh2,K,rjac,omega=1.0;
	H=1.0;
	h=1.0/(n-1); 
	h2i=1.0/(h*h);
	foh2=-4.0*h2i;
	K=1.0;
	rjac=/*0.99555*/0.9999;
	for (jpost=1;jpost<=NPOST;jpost++)
	{ jsw=1;
	for (ipass=1;ipass<=2;ipass++)
	{
		lsw=jsw;
		for(i=2;i<n;i++)
		{
			for (j=lsw+1;j<n;j+=2)
			{
				        v1=u[i][j+1]-u[i][j];
						v2=u[i][j]-u[i][j-1];
						v4=u[i+1][j]-u[i][j];
						v5=u[i][j]-u[i-1][j];
						v6=u[i+1][j+1]-u[i-1][j+1];
						v7=u[i-1][j-1]-u[i+1][j-1];
						v3=-4.0*pow(1.0*(i-1)*q[j]/H,2)-4.0*(pow(1.0*(i-1),2) -pow(1.0*(i-1)*(q[j+1]-q[j-1])/H,2)*pow(1.0*(i-1),2))-2.0*cos(2.0*u[i][j]);
						v8=q[j+1]-q[j];
						v9=q[j]-q[j-1];
						res=2.0*pow(1.0*(i-1)*q[j]/H,2)*(v1-v2)+(2.0*pow(1.0*(i-1),2)+0.5*pow(1.0*(i-1)*(q[j+1]-q[j-1])/H,2)*pow(1.0*(i-1),2))*(v4-v5)+(1.0*(i-1)-pow(1.0*(i-1),3)*q[j]*(q[j+1]+q[j-1]-2.0*q[j])/pow(H,2)+0.5*pow(1.0*(i-1),3)*pow(q[j+1]-q[j-1],2)/pow(H,2))*(v4+v5)-(0.5*pow(1.0*(i-1),3)*q[j]*(q[j+1]-q[j-1])/pow(H,2))*(v6+v7)-sin(2.0*u[i][j]);
						u[i][j] -= omega*res/v3;
							for (j=2;j<=1+b*(int)(R0/H);j++)
	{
		u[a+b][j]=0.5*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))+0.5*Pi;
	}
	for (i=a+1;i<=a+b;i++)
	{
		u[i][1+b*(int)(R0/H)]=0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))+0.5*Pi;
	}
	u[a][1+b*(int)(R0/H)]=0.75*Pi;
	for (i=a-b;i<=a-1;i++)
    {
		u[i][1+b*(int)(R0/H)]=0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))+Pi;
	}
	for (j=2;j<=1+b*(int)(R0/H);j++)
    {
		u[a-b][j]=0.5*atan(H*(y[j]-y[1])/(R0*(x[a-b]-x[a])))+Pi;
	}
	for (i=n-b;i<=n-1;i++)
	{
		u[i][n-b*(int)(R0/H)]=Pi-atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])));
	}
	for (j=n-b*(int)(R0/H);j<=n-1;j++)
    {
		u[n-b][j]=Pi-atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])));
	}
			}
			lsw=3-lsw;
		}
		jsw=3-jsw;
		omega=(jpost==1&&ipass==1?1.0/(1.0-0.5*rjac*rjac):1.0/(1.0-0.25*rjac*rjac*omega));
	}
	}

	
}