#include <stdio.h>
#include <math.h>
#define Pi 3.141592653589793
#define Up 1.570796326794897
#define Down1 1.570796326794897
#define Down2 1.570796326794897
#define Right 1.570796326794897
#define Left 1.570796326794897
#define S 0.00000000000001
void relax2(double **u, double **rhs, int n, double m,int a,double *q) 
{ 
	int i,ipass,isw,j,jsw=1; 
	double h,res,v1,v2,v3,v4,v5,v6,v7,v8,v9,H,h2i,foh2,K;
	H=1.0;
	h=1.0/(n-1); 
	h2i=1.0/(h*h);
	foh2=-4.0*h2i;
	K=1.0;
	/*for (i=1;i<=n;i++)
		for (j=1;j<=n;j++)
		{
			if (u[i][j]<0.0)
				u[i][j]=0.0;
			
			
			if (u[i][j]>Pi)
				u[i][j]=Pi;
			
		}*/
	for (i=1;i<=a;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=n;i++)
		u[i][1]=Down2;
	for (i=1;i<=n;i++)
		u[i][n]=Up-0.5*m*Pi;

	for (j=1;j<=n;j++)
		{
			u[n][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	}

	for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) 
	{ 
		isw=jsw; 
		for (j=2;j<=n-1;j++,isw=3-isw) 
		{ 
			for (i=isw+1;i<n;i+=2) 
			{ 
	

				
						v1=u[i][j+1]-u[i][j];
						v2=u[i][j]-u[i][j-1];
						v4=u[i+1][j]-u[i][j];
						v5=u[i][j]-u[i-1][j];
						v6=u[i+1][j+1]-u[i-1][j+1];
						v7=u[i-1][j-1]-u[i+1][j-1];
						/*v3=-4.0*pow(1.0*(i-1)*q[j]/H,2)-4.0*(pow(1.0*(i-1),2) -pow(1.0*(i-1)*(q[j+1]-q[j-1])/H,2)*pow(1.0*(i-1),2))-2.0*cos(2.0*u[i][j]);*/
						v8=q[j+1]-q[j];
						v9=q[j]-q[j-1];

				
				res=/*2.0*pow(1.0*(i-1)*q[j]/H,2)*(v1-v2)+(2.0*pow(1.0*(i-1),2)+0.5*pow(1.0*(i-1)*(q[j+1]-q[j-1])/H,2)*pow(1.0*(i-1),2))*(v4-v5)+(1.0*(i-1)-pow(1.0*(i-1),3)*q[j]*(q[j+1]+q[j-1]-2.0*q[j])/pow(H,2)+0.5*pow(1.0*(i-1),3)*pow(q[j+1]-q[j-1],2)/pow(H,2))*(v4+v5)-(0.5*pow(1.0*(i-1),3)*q[j]*(q[j+1]-q[j-1])/pow(H,2))*(v6+v7)-sin(2.0*u[i][j])-rhs[i][j];*/
				
				K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*q[j]*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2)-sin(2.0*u[i][j])*0.25*pow(1.0*(i-1),4)*pow(v8+v9,2)-cos(2.0*u[i][j])*pow(1.0*(i-1),3)*(v8+v9))*0.25*pow(v4+v5,2)
				+(K-1.0)*(0.5*(v8+v9)*q[j]*sin(2.0*u[i][j])*pow(1.0*(i-1),3)+q[j]*cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)+pow(1.0*(i-1),3)*(0.5*pow(v8+v9,2)-(v8-v9)*q[j]))-(K-1.0)*pow(1.0*(i-1),2)*0.75*(v8+v9))*(v4+v5)
				+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j]*0.5*(v1+v2)-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),4)*0.25*pow(v8+v9,2)+(K-1)*pow(1.0*(i-1),3)*0.5*(v8+v9)*sin(2.0*u[i][j]))*(v4-v5)
				-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1)*q[j],2)*(v1-v2)+((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*(v8+v9)*q[j]+(K-1.0)*q[j]*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7)-rhs[i][j];
				
				v3=2.0*K*cos(2.0*u[i][j])-0.5*(K-1.0)*pow(q[j]*(i-1)*(v1+v2),2)*cos(2.0*u[i][j])+(K-1.0)*(2.0*cos(2.0*u[i][j])*pow(1.0*(i-1),2)-0.5*pow(1.0*(i-1),4)*pow(v8+v9,2)*cos(2.0*u[i][j])+2.0*pow(1.0*(i-1),3)*(v8+v9)*sin(2.0*u[i][j]))*0.25*pow(v4+v5,2)
					+(K-1.0)*(pow(1.0*(i-1),3)*(v8+v9)*q[j]*cos(2.0*u[i][j])-2.0*(i-1)*(i-1)*q[j]*sin(2.0*u[i][j]))*0.5*(v1+v2)*(v4+v5)+(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*(1.0*(i-1)+pow(1.0*(i-1),3)*(0.5*pow(v8+v9,2)-(v8-v9)*q[j]))*(v4+v5)
					+(K-1.0)*cos(2.0*u[i][j])*(i-1)*q[j]*(v1+v2)-2.0*(-(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),2)+(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),4)*0.25*pow(v8+v9,2)+(K-1.0)*pow(1.0*(i-1),3)*(v8+v9)*cos(2.0*u[i][j]))*(v4-v5)
					+4.0*((K*pow(cos(1.0*u[i][j]),2)+1.0*pow(sin(1.0*u[i][j]),2))*pow(1.0*(i-1),2)+(K*pow(sin(1.0*u[i][j]),2)+1.0*pow(cos(1.0*u[i][j]),2))*pow(1.0*(i-1),4)*0.25*pow(v8+v9,2)+(K-1.0)*pow(1.0*(i-1),3)*0.5*(v8+v9)*sin(2.0*u[i][j]))
					-2.0*(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1)*q[j],2)*(v1-v2)+4.0*(K*pow(sin(1.0*u[i][j]),2)+1.0*pow(cos(1.0*u[i][j]),2))*pow(1.0*(i-1)*q[j],2)+((K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),3)*(v8+v9)*q[j]+2.0*(K-1.0)*q[j]*pow(1.0*(i-1),2)*cos(2.0*u[i][j]))*0.5*(v6+v7);
				u[i][j] -= res/v3; 
			}
}
	}
}