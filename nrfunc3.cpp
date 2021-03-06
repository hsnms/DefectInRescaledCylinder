#include <math.h>
#include <stdio.h>
#include "Header2.h"
#include "Header3.h"
/*#define EPI 1.0e-7*/ 
#define Pi 3.141592653589793
void **nrfunc3(double **f, double **u, int n, double *q, int a,int b,double R0, double K)
{
	double h,v1,v2,v4,v5,v6,v7,v8,v9,H;
	int i,j;
	h=1.0/(n-1); 	
	H=1.0;
	for(j=1;j<=n;j++)
		for(i=1;i<=n;i++)
			f[i][j]=0.0;


		for(i=2;i<n-b;i++)
	{
		for (j=b*(int)(R0/H)+2;j<n;j++)
		{
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];	        
			/*f[i][j]=pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
				
		}
	}

	
		for(i=2;i<a-b;i++)
	{
		for (j=2;j<=1+b*(int)(R0/H);j++)
		{
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];	        
			/*f[i][j]=pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
				
		}
	}	

		
		for(i=a+b+1;i<n;i++)
	{
		for (j=2;j<=1+b*(int)(R0/H);j++)
		{
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];	        
			/*f[i][j]=pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
				
		}
	}

		
		for(i=n-b;i<n;i++)
	{
		for (j=b*(int)(R0/H)+2;j<n-b*(int)(R0/H);j++)
		{
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];	        
			/*f[i][j]=pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
				
		}
	}

	for (i=2;i<=a-b-1;i++)
	{
			v1=u[i][2]-u[i][1];
		    v2=u[i][2]-u[i][1];
			v4=u[i+1][1]-u[i][1];
			v5=u[i][1]-u[i-1][1];
			v8=q[2]-q[1];
			v9=q[2]-q[1];	
			j=1;
		/*f[i][1]=pow(sin(u[i][1]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][1])*(v4+v5)/(2.0*h)-pow(sin(u[i][1]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			

	}

		for (i=a+b+1;i<=n-1;i++)
	{
			v1=u[i][2]-u[i][1];
		    v2=u[i][2]-u[i][1];
			v4=u[i+1][1]-u[i][1];
			v5=u[i][1]-u[i-1][1];
			v8=q[2]-q[1];
			v9=q[2]-q[1];	
			j=1;
		/*f[i][1]=pow(sin(u[i][1]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][1])*(v4+v5)/(2.0*h)-pow(sin(u[i][1]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			

	}

	for (i=2;i<=n-b-1;i++)
	{
			v1=u[i][n]-u[i][n-1];
		    v2=u[i][n]-u[i][n-1];
			v4=u[i+1][n]-u[i][n];
			v5=u[i][n]-u[i-1][n];
			v8=q[n]-q[n-1];
			v9=q[n]-q[n-1];	  
			j=n;
		/*f[i][n]=pow(sin(u[i][n]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][n])*(v4+v5)/(2.0*h)-pow(sin(u[i][n]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			

	}
	for (j=2;j<=n-1;j++)
	{
		    v1=u[1][j+1]-u[1][j];
		    v2=u[1][j]-u[1][j-1];
			v4=u[2][j]-u[1][j];
			v5=u[2][j]-u[1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];		       
			i=1;
		/*f[1][j]=sin(2.0*u[1][j])*(v4+v5)/(2.0*h)-pow(sin(u[1][j]),2)*q[j]*(v1+v2)/(h);*/
			f[i][j]=((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			
	}

	for (j=2;j<=n-b*(int)(R0/H)-1;j++)
	{
		    v1=u[n][j+1]-u[n][j];
		    v2=u[n][j]-u[n][j-1];
			v4=u[n][j]-u[n-1][j];
			v5=u[n][j]-u[n-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];		        
			i=n;
		/*f[n][j]=pow(sin(u[n][j]),2)/(1.0*(n-1)*h)+sin(2.0*u[n][j])*(v4+v5)/(2.0*h)-pow(sin(u[n][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(n-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(n-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			
	}

		for (j=2;j<=b*(int)(R0/H);j++)
	{
		    v1=u[a-b][j+1]-u[a-b][j];
		    v2=u[a-b][j]-u[a-b][j-1];
			v4=u[a-b][j]-u[a-b-1][j];
			v5=u[a-b][j]-u[a-b-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];		        
			i=a-b;
		/*f[i][j]=pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
	}

		for (j=n-b*(int)(R0/H)+1;j<=n-1;j++)
	{
		    v1=u[n-b][j+1]-u[n-b][j];
		    v2=u[n-b][j]-u[n-b][j-1];
			v4=u[n-b][j]-u[n-b-1][j];
			v5=u[n-b][j]-u[n-b-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];		        
			i=n-b;
		/*f[i][j]=pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
	}

	for (j=2;j<=b*(int)(R0/H);j++)
	{
		    v1=u[a+b][j+1]-u[a+b][j];
		    v2=u[a+b][j]-u[a+b][j-1];
			v4=u[a+b+1][j]-u[a+b][j];
			v5=u[a+b+1][j]-u[a+b][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];		        
			i=a+b;
		/*f[i][j]=pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
	}


		for (i=a-b+1;i<=a+b-1;i++)
	{
		    v1=u[i][b*(int)(R0/H)+2]-u[i][b*(int)(R0/H)+1];
		    v2=u[i][b*(int)(R0/H)+2]-u[i][b*(int)(R0/H)+1];
			v4=u[i+1][b*(int)(R0/H)+1]-u[i][b*(int)(R0/H)+1];
			v5=u[i][b*(int)(R0/H)+1]-u[i-1][b*(int)(R0/H)+1];
			v8=q[b*(int)(R0/H)+2]-q[b*(int)(R0/H)+1];
			v9=q[b*(int)(R0/H)+2]-q[b*(int)(R0/H)+1];		        
			j=b*(int)(R0/H)+1;
		/*f[i][j]=pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
	}

		for (i=n-b+1;i<=n-1;i++)
	{
		    v1=u[i][n-b*(int)(R0/H)]-u[i][n-b*(int)(R0/H)-1];
		    v2=u[i][n-b*(int)(R0/H)]-u[i][n-b*(int)(R0/H)-1];
			v4=u[i+1][n-b*(int)(R0/H)]-u[i][n-b*(int)(R0/H)];
			v5=u[i][n-b*(int)(R0/H)]-u[i-1][n-b*(int)(R0/H)];
			v8=q[n-b*(int)(R0/H)]-q[n-b*(int)(R0/H)-1];
			v9=q[n-b*(int)(R0/H)]-q[n-b*(int)(R0/H)-1];	        
			j=n-b*(int)(R0/H);
		/*f[i][j]=pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
	}
            v1=u[1][2]-u[1][1];
		    v2=u[1][2]-u[1][1];
			v4=u[2][1]-u[1][1];
			v5=u[2][1]-u[1][1];
			v8=q[2]-q[1];
			v9=q[2]-q[1];		       
			i=1;j=1;
		/*f[1][j]=sin(2.0*u[1][j])*(v4+v5)/(2.0*h)-pow(sin(u[1][j]),2)*q[j]*(v1+v2)/(h);*/
			f[i][j]=((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
	  
			v1=u[1][n]-u[1][n-1];
		    v2=u[1][n]-u[1][n-1];
			v4=u[2][n]-u[1][n];
			v5=u[2][n]-u[1][n];
			v8=q[n]-q[n-1];
			v9=q[n]-q[n-1];		       
			i=1;j=n;
		/*f[1][j]=sin(2.0*u[1][j])*(v4+v5)/(2.0*h)-pow(sin(u[1][j]),2)*q[j]*(v1+v2)/(h);*/
			f[i][j]=((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			
	        v1=u[n][2]-u[n][1];
		    v2=u[n][2]-u[n][1];
			v4=u[n][1]-u[n-1][1];
			v5=u[n][1]-u[n-1][1];
			v8=q[2]-q[1];
			v9=q[2]-q[1];		        
			i=n;j=1;
		/*f[n][j]=pow(sin(u[n][j]),2)/(1.0*(n-1)*h)+sin(2.0*u[n][j])*(v4+v5)/(2.0*h)-pow(sin(u[n][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(n-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(n-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2);*/
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
	
			/*strange, shouldn't have this part*/			/* v1=u[n][n]-u[n][n-1];
		    v2=u[n][n]-u[n][n-1];
			v4=u[n][n]-u[n-1][n];
			v5=u[n][n]-u[n-1][n];
			v8=q[n]-q[n-1];
			v9=q[n]-q[n-1];		        
			i=n;j=n;
	
			f[i][j]=K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);*/
	
            
	        v1=u[a-b][2]-u[a-b][1];
		    v2=u[a-b][2]-u[a-b][1];
			v4=u[a-b][1]-u[a-b-1][1];
			v5=u[a-b][1]-u[a-b-1][1];
			v8=q[2]-q[1];
			v9=q[2]-q[1];
			i=a-b;j=1;
			f[a-b][1]=/*pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2)*/
				K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
	
           
			v1=u[a+b][2]-u[a+b][1];
		    v2=u[a+b][2]-u[a+b][1];
			v4=u[a+b+1][1]-u[a+b][1];
			v5=u[a+b+1][1]-u[a+b][1];
			v8=q[2]-q[1];
			v9=q[2]-q[1];	
			i=a+b;j=1;
			f[i][j]=/*pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2)*/
				K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			
			
			v1=u[n-b][n]-u[n-b][n-1];
		    v2=u[n-b][n]-u[n-b][n-1];
			v4=u[n-b][n]-u[n-b-1][n];
			v5=u[n-b][n]-u[n-b-1][n];
			v8=q[n]-q[n-1];
			v9=q[n]-q[n-1];	 
			i=n-b;j=n;
			f[i][j]=/*pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2)*/
				K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			
			
			v1=u[n][n-b*(int)(R0/H)]-u[n][n-b*(int)(R0/H)-1];
		    v2=u[n][n-b*(int)(R0/H)]-u[n][n-b*(int)(R0/H)-1];
			v4=u[n][n-b*(int)(R0/H)]-u[n-1][n-b*(int)(R0/H)];
			v5=u[n][n-b*(int)(R0/H)]-u[n-1][n-b*(int)(R0/H)];
			v8=q[n-b*(int)(R0/H)]-q[n-b*(int)(R0/H)-1];
			v9=q[n-b*(int)(R0/H)]-q[n-b*(int)(R0/H)-1];
			i=n;j=n-b*(int)(R0/H);
			f[i][j]=/*pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2)*/
				K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			
			i=a-b;j=b*(int)(R0/H)+1;
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];	  
			f[i][j]=/*pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2)*/
				K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			
			i=a+b;j=b*(int)(R0/H)+1;
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];	  
			f[i][j]=/*pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2)*/
				K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			
			i=n-b;j=n-b*(int)(R0/H);
			v1=u[i][j+1]-u[i][j];
		    v2=u[i][j]-u[i][j-1];
			v4=u[i+1][j]-u[i][j];
			v5=u[i][j]-u[i-1][j];
			v8=q[j+1]-q[j];
			v9=q[j]-q[j-1];	  
			f[i][j]=/*pow(sin(u[i][j]),2)/(1.0*(i-1)*h)+sin(2.0*u[i][j])*(v4+v5)/(2.0*h)-pow(sin(u[i][j]),2)*q[j]*(v1+v2)/(h)+(1.0*(i-1)*h)*pow((v4+v5)/(2.0*h),2)+(1.0*(i-1)*h)*pow(q[j],2)*pow((v1+v2)/(2.0*h),2)*/
				K*pow(sin(u[i][j]),2)/(h*(i-1))+((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(i-1)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*0.25*pow(v8+v9,2)+(K-1.0)*sin(2.0*u[i][j])*pow(1.0*(i-1),2)*0.5*(v8+v9))*pow(v4+v5,2)/(4.0*h)
				+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*(i-1)*pow(q[j]*(v1+v2),2)/(4.0*h)-((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),2)*(v8+v9)*q[j]+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j])*(v1+v2)*(v4+v5)/(4.0*h)
				+(K*sin(2.0*u[i][j])-K*(cos(2.0*u[i][j])-1.0)*(i-1)*0.5*(v8+v9))*(v4+v5)/(2.0*h)+K*(cos(2.0*u[i][j])-1.0)*q[j]*(v1+v2)/(2.0*h);
			
	return 0;
}
