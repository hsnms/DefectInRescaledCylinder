#include <math.h>
#include <stdio.h>
#define Pi 3.141592653589793
#define S 0.00000000000001
double Stat3 (double **u,int n,double *q,int b, double R0, double K)
{
	double zero,v1,v2,v4,v5,v6,v7,s1,s2,H,h2i,h,v8,v9;
	int i,j,k1,k2;
	zero=0.0;h=1.0/(n-1); 
	H=1.0;h2i=1.0/(h*h);
	for(i=2;i<n;i++)
	{
		s2=0.0;
		for (j=b*(int)(R0/H)+2;j<n-b*(int)(R0/H);j++)
		{
			v1=(u[i][j+1]-u[i][j]);
			v2=(u[i][j]-u[i][j-1]);
			v4=(u[i+1][j]-u[i][j]);
			v5=(u[i][j]-u[i-1][j]);
			v6=u[i+1][j+1]-u[i-1][j+1];
			v7=u[i-1][j-1]-u[i+1][j-1];  
			v8=q[j+1]-q[j];
		    v9=q[j]-q[j-1];           
								 
				/*s1=2.0*pow(1.0*(i-1)*q[j]/H,2)*(v1-v2)+(2.0*pow(1.0*(i-1),2)+0.5*pow(1.0*(i-1)*(q[j+1]-q[j-1])/H,2)*pow(1.0*(i-1),2))*(v4-v5)+(1.0*(i-1)-pow(1.0*(i-1),3)*q[j]*(q[j+1]+q[j-1]-2.0*q[j])/pow(H,2)+0.5*pow(1.0*(i-1),3)*pow(q[j+1]-q[j-1],2)/pow(H,2))*(v4+v5)-(0.5*pow(1.0*(i-1),3)*q[j]*(q[j+1]-q[j-1])/pow(H,2))*(v6+v7)-sin(2.0*u[i][j]);*/
			s1=	K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*q[j]*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2)-sin(2.0*u[i][j])*0.25*pow(1.0*(i-1),4)*pow(v8+v9,2)-cos(2.0*u[i][j])*pow(1.0*(i-1),3)*(v8+v9))*0.25*pow(v4+v5,2)
				+(K-1.0)*(0.5*(v8+v9)*q[j]*sin(2.0*u[i][j])*pow(1.0*(i-1),3)+q[j]*cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)+pow(1.0*(i-1),3)*(0.5*pow(v8+v9,2)-(v8-v9)*q[j]))-(K-1.0)*pow(1.0*(i-1),2)*0.75*(v8+v9))*(v4+v5)
				+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j]*0.5*(v1+v2)-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),4)*0.25*pow(v8+v9,2)+(K-1)*pow(1.0*(i-1),3)*0.5*(v8+v9)*sin(2.0*u[i][j]))*(v4-v5)
				-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1)*q[j],2)*(v1-v2)+((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*(v8+v9)*q[j]+(K-1.0)*q[j]*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
			s2=s2+fabs(s1);
			
		}
		zero=zero+s2;
	}
		for(i=b+2;i<n;i++)
	{
		s2=0.0;
		for (j=2;j<=1+b*(int)(R0/H);j++)
		{
			v1=(u[i][j+1]-u[i][j]);
			v2=(u[i][j]-u[i][j-1]);
			v4=(u[i+1][j]-u[i][j]);
			v5=(u[i][j]-u[i-1][j]);
			v6=u[i+1][j+1]-u[i-1][j+1];
			v7=u[i-1][j-1]-u[i+1][j-1];  
			v8=q[j+1]-q[j];
		    v9=q[j]-q[j-1];           
								 
				/*s1=2.0*pow(1.0*(i-1)*q[j]/H,2)*(v1-v2)+(2.0*pow(1.0*(i-1),2)+0.5*pow(1.0*(i-1)*(q[j+1]-q[j-1])/H,2)*pow(1.0*(i-1),2))*(v4-v5)+(1.0*(i-1)-pow(1.0*(i-1),3)*q[j]*(q[j+1]+q[j-1]-2.0*q[j])/pow(H,2)+0.5*pow(1.0*(i-1),3)*pow(q[j+1]-q[j-1],2)/pow(H,2))*(v4+v5)-(0.5*pow(1.0*(i-1),3)*q[j]*(q[j+1]-q[j-1])/pow(H,2))*(v6+v7)-sin(2.0*u[i][j]);*/
			s1=	K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*q[j]*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2)-sin(2.0*u[i][j])*0.25*pow(1.0*(i-1),4)*pow(v8+v9,2)-cos(2.0*u[i][j])*pow(1.0*(i-1),3)*(v8+v9))*0.25*pow(v4+v5,2)
				+(K-1.0)*(0.5*(v8+v9)*q[j]*sin(2.0*u[i][j])*pow(1.0*(i-1),3)+q[j]*cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)+pow(1.0*(i-1),3)*(0.5*pow(v8+v9,2)-(v8-v9)*q[j]))-(K-1.0)*pow(1.0*(i-1),2)*0.75*(v8+v9))*(v4+v5)
				+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j]*0.5*(v1+v2)-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),4)*0.25*pow(v8+v9,2)+(K-1)*pow(1.0*(i-1),3)*0.5*(v8+v9)*sin(2.0*u[i][j]))*(v4-v5)
				-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1)*q[j],2)*(v1-v2)+((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*(v8+v9)*q[j]+(K-1.0)*q[j]*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
			s2=s2+fabs(s1);
			
		}
		zero=zero+s2;
	}
	for(i=2;i<n-b;i++)
	{
		s2=0.0;
		for (j=n-b*(int)(R0/H);j<n;j++)
		{
			v1=(u[i][j+1]-u[i][j]);
			v2=(u[i][j]-u[i][j-1]);
			v4=(u[i+1][j]-u[i][j]);
			v5=(u[i][j]-u[i-1][j]);
			v6=u[i+1][j+1]-u[i-1][j+1];
			v7=u[i-1][j-1]-u[i+1][j-1];  
			v8=q[j+1]-q[j];
		    v9=q[j]-q[j-1];           
								 
				/*s1=2.0*pow(1.0*(i-1)*q[j]/H,2)*(v1-v2)+(2.0*pow(1.0*(i-1),2)+0.5*pow(1.0*(i-1)*(q[j+1]-q[j-1])/H,2)*pow(1.0*(i-1),2))*(v4-v5)+(1.0*(i-1)-pow(1.0*(i-1),3)*q[j]*(q[j+1]+q[j-1]-2.0*q[j])/pow(H,2)+0.5*pow(1.0*(i-1),3)*pow(q[j+1]-q[j-1],2)/pow(H,2))*(v4+v5)-(0.5*pow(1.0*(i-1),3)*q[j]*(q[j+1]-q[j-1])/pow(H,2))*(v6+v7)-sin(2.0*u[i][j]);*/
			s1=	K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*q[j]*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2)-sin(2.0*u[i][j])*0.25*pow(1.0*(i-1),4)*pow(v8+v9,2)-cos(2.0*u[i][j])*pow(1.0*(i-1),3)*(v8+v9))*0.25*pow(v4+v5,2)
				+(K-1.0)*(0.5*(v8+v9)*q[j]*sin(2.0*u[i][j])*pow(1.0*(i-1),3)+q[j]*cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)+pow(1.0*(i-1),3)*(0.5*pow(v8+v9,2)-(v8-v9)*q[j]))-(K-1.0)*pow(1.0*(i-1),2)*0.75*(v8+v9))*(v4+v5)
				+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j]*0.5*(v1+v2)-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),4)*0.25*pow(v8+v9,2)+(K-1)*pow(1.0*(i-1),3)*0.5*(v8+v9)*sin(2.0*u[i][j]))*(v4-v5)
				-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1)*q[j],2)*(v1-v2)+((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*(v8+v9)*q[j]+(K-1.0)*q[j]*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7);
			s2=s2+fabs(s1);
			
		}
		zero=zero+s2;
	}



	
	zero=zero/((n-2)*(n-2)-2*(b-1)*(b*(int)(R0/H)-1));
	return zero;
}