#include <math.h> 
#include <stdio.h>
#include "Header1.h"
#define MAXIT 100 
#define RS 0.1
#define S 0.00000000000001
#define Pi 3.141592653589793
#define A pow(sin(u[n][j]),2)+v4*sin(2.0*u[n][j])/h-pow(sin(u[n][j]),2)*(v1+v2)/(h*H)
#define B v4*pow(sin(u[n][j])/h,2)*2.0/H
#define C pow(1.0/h,2)*pow(v4,2)+pow((v1+v2)/(2.0*h*H),2)
#define D pow(v4/(h*H),2)/pow(h,2)
#define E -(v1+v2)*v4/(pow(h,3)*pow(H,2))
#define F Pi*s*(2.0*pow((1+((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1])/(h*H)),0.5)-((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1])*pow((1+((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1])/(h*H)),-0.5)/(h*H)+0.5*(q[j+1]+q[j-1]-2.0*q[j])*q[j]*pow((1+((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1])/(h*H)),-1.5)/pow(h*H,2))
double rtsafe( double H, double h,int n,double x1, double x2, double xcc, double rrhs, double **u, double *q) 

{ 
	double df,dx,dxold,f,fh,fl,max,min,x3,xmin,fm,xmax; 
	double temp,xh,xl,rts;
	int i;
	void fanddf (double x,double *f,double *df,double **u, double *q, double H, double h,int n, double rrhs);
fanddf(x1, &fl, &df,u,q,H,h,n,rrhs);
fanddf(x2, &fh, &df,u,q,H,h,n,rrhs);
/*printf("fl=%f\tfh=%f\n",fl,fh);*/
if (fl > 0.0 && fh > 0.0) 
	{ 
		if (fl>=fh) {max=fl; min=fh;xmin=x2;}
		else {max=fh;min=fl;xmin=x1;}
		for (x3=x1;x3<=x2;x3=x3+RS)
		{
			
			fanddf(x3, &fm, &df,u,q,H,h,n,rrhs);
			if (fm<=min)
			{
				min=fm;
				xmin=x3;
			}
		} return xmin;
/*printf("---min=%f\txmin=%f---\n",min,xmin);*/
		/*printf("fl=%f\tfh=%f\n",fl,fh);
nrerror2("Root must be bracketed in rtsafe");*/
}
if (fl < 0.0 && fh < 0.0)
{  
	if (fl>=fh) {max=fl; min=fh;xmax=x1;}
		else {max=fh;min=fl;xmax=x2;}
		for (x3=x1;x3<=x2;x3=x3+RS)
		{
			
			fanddf(x3, &fm, &df,u,q,H,h,n,rrhs);
			if (fm>=max)
			{
				max=fm;
				xmax=x3;
			}
		} return xmax;
/*printf("---max=%f\txmax=%f---\n",max,xmax);*/
	
		/*nrerror2("Root must be bracketed in rtsafe");*/
}

if (fl == 0.0) return x1; 
if (fh == 0.0) return x2; 
if (fl < 0.0)
{
	xl=x1; xh=x2; 
} 
else 
{ 
	xh=x1; 
	xl=x2; 
} 
rts=0.5*(x1+x2); 
dxold=fabs(x2-x1); 
dx=dxold; 
fanddf(rts, &f, &df,u,q,H,h,n,rrhs);
for (i=1;i<=MAXIT;i++) 
{
	if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)  || (fabs(2.0*f) > fabs(dxold*df))) 
	{ 
		dxold=dx; 
		dx=0.5*(xh-xl); 
		rts=xl+dx; 
		if (xl == rts) return rts;
	}
	else 
	{ 
		dxold=dx;
		dx=f/df; 
		temp=rts; 
		rts -= dx; 
		if (temp == rts) return rts; 
	} 
	if (fabs(dx) < xcc) return rts; 

	if (f < 0.0) 
		xl=rts;
	else xh=rts;
} 
nrerror2("Maximum number of iterations exceeded in rtsafe"); 
return 0.0; 
}