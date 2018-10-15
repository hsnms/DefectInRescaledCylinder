#include <math.h>
#include <stdio.h>
#include "Header2.h"

void slvsm2(int n, double **u, double **rhs, double m, int a, double *q) 
{ 
	double H,h,x1,x2,xcc,rrhs,fact,disc;
	double rtsafe( double H, double h,int n,double x1, double x2,double xcc, double rrhs,double **u, double *q) ;
	H=1.0;
	h=/*1.0/(n-1)*/0.5; 
	fact=2.0/(h*h);
	disc=sqrt(fact*fact+rhs[2][2]);
	fill02(u,n,m,a); 
	
	xcc=0.001;
	/*u[2][2]=-rhs[2][2]*pow(H/(2.0*q[2]),2);*/
	/*u[2][2]=*//*-rhs[2][2]/(fact+disc)*//*(2.0-sqrt(4.0+pow(h,4)*rhs[2][2]*(u[3][2]+u[1][2]+u[2][3]+u[2][1])))/pow(h,2);*/
	x1=u[1][2];
	x2=u[3][2];
	rrhs=rhs[2][2];
	u[2][2]=rtsafe( H, h,n, x1,x2, xcc,rrhs,u,q);
	/*printf("u[2][2]=%f\n",u[2][2]);*/
	
}