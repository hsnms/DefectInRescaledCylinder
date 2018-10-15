#include <math.h> 
#define S 0.00000000000001
#define Pi 3.141592653589793
#define A pow(sin(u[n][j]),2)+v4*sin(2.0*u[n][j])/h-pow(sin(u[n][j]),2)*(v1+v2)/(h*H)
#define B v4*pow(sin(u[n][j])/h,2)*2.0/H
#define C pow(1.0/h,2)*pow(v4,2)+pow((v1+v2)/(2.0*h*H),2)
#define D pow(v4/(h*H),2)/pow(h,2)
#define E -(v1+v2)*v4/(pow(h,3)*pow(H,2))
#define F 2.0*Pi*s*H*pow((1+pow(((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1])/(h*H),2)),-1.5)*(1+pow(((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1])/(h*H),2)-q[j]*(q[j+1]+q[j-1]-2.0*q[j])/pow(h*H,2))
void fanddf (double x, double *f,double *df,double **u, double *q, double H, double h,int n, double rrhs)
{
	double v1,v2,v4,v5,v6,v7,K,v8,v9;
	int i,j;
	K=1.0;
	                    v1=u[2][3]-x;
						v2=x-u[2][1];
						v4=u[3][2]-x;
						v5=x-u[1][2];
						v6=u[3][3]-u[1][3];
						v7=u[1][1]-u[3][1];
    	                v8=q[3]-q[3];
						v9=q[3]-q[1];
		i=2;j=2;u[2][2]=x;				
		*f=2.0*pow(1.0*q[2]/H,2)*(v1-v2)+(2.0*pow(1.0,2)+0.5*pow(1.0*(q[3]-q[1])/H,2)*pow(1.0,2))*(v4-v5)+(1.0-pow(1.0,3)*q[2]*(q[3]+q[1]-2.0*q[2])/pow(H,2)+0.5*pow(1.0,3)*pow(q[3]-q[1],2)/pow(H,2))*(v4+v5)-(0.5*pow(1.0,3)*q[2]*(q[3]-q[1])/pow(H,2))*(v6+v7)-sin(2.0*u[2][2])-rrhs;
	/*K*sin(2.0*u[i][j])-(K-1.0)*0.25*sin(2.0*u[i][j])*pow((i-1)*q[j]*(v1+v2),2)+(K-1.0)*(sin(2.0*u[i][j])*pow(1.0*(i-1),2)-sin(2.0*u[i][j])*0.25*pow(1.0*(i-1),4)*pow(v8+v9,2)-cos(2.0*u[i][j])*pow(1.0*(i-1),3)*(v8+v9))*0.25*pow(v4+v5,2)
				+(K-1.0)*(0.5*(v8+v9)*q[j]*sin(2.0*u[i][j])*pow(1.0*(i-1),3)+q[j]*cos(2.0*u[i][j])*pow(1.0*(i-1),2))*0.5*(v4+v5)*(v1+v2)+(-(K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*(1.0*(i-1)+pow(1.0*(i-1),3)*(0.5*pow(v8+v9,2)-(v8-v9)*q[j]))-(K-1.0)*pow(1.0*(i-1),2)*0.75*(v8+v9))*(v4+v5)
				+(K-1.0)*sin(2.0*u[i][j])*(i-1)*q[j]*0.5*(v1+v2)-2.0*((K*pow(cos(u[i][j]),2)+1.0*pow(sin(u[i][j]),2))*pow(1.0*(i-1),2)+(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),4)*0.25*pow(v8+v9,2)+(K-1)*pow(1.0*(i-1),3)*0.5*(v8+v9)*sin(2.0*u[i][j]))*(v4-v5)
				-2.0*(K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1)*q[j],2)*(v1-v2)+((K*pow(sin(u[i][j]),2)+1.0*pow(cos(u[i][j]),2))*pow(1.0*(i-1),3)*(v8+v9)*q[j]+(K-1.0)*q[j]*pow(1.0*(i-1),2)*sin(2.0*u[i][j]))*0.5*(v6+v7)-rrhs;*/
		*df=-4.0*pow(1.0*q[2]/H,2)-4.0*(pow(1.0,2) -pow(1.0*(q[3]-q[1])/H,2)*pow(1.0,2))-2.0*cos(2.0*u[2][2]);
			/*2.0*K*cos(2.0*u[i][j])-0.5*(K-1.0)*pow(q[j]*(i-1)*(v1+v2),2)*cos(2.0*u[i][j])+(K-1.0)*(2.0*cos(2.0*u[i][j])*pow(1.0*(i-1),2)-0.5*pow(1.0*(i-1),4)*pow(v8+v9,2)*cos(2.0*u[i][j])+2.0*pow(1.0*(i-1),3)*(v8+v9)*sin(2.0*u[i][j]))*0.25*pow(v4+v5,2)
					+(K-1.0)*(pow(1.0*(i-1),3)*(v8+v9)*q[j]*cos(2.0*u[i][j])-2.0*(i-1)*(i-1)*q[j]*sin(2.0*u[i][j]))*0.5*(v1+v2)*(v4+v5)+(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*(1.0*(i-1)+pow(1.0*(i-1),3)*(0.5*pow(v8+v9,2)-(v8-v9)*q[j]))*(v4+v5)
					+(K-1.0)*cos(2.0*u[i][j])*(i-1)*q[j]*(v1+v2)-2.0*(-(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),2)+(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),4)*0.25*pow(v8+v9,2)+(K-1.0)*pow(1.0*(i-1),3)*(v8+v9)*cos(2.0*u[i][j]))*(v4-v5)
					+4.0*((K*pow(cos(1.0*u[i][j]),2)+1.0*pow(sin(1.0*u[i][j]),2))*pow(1.0*(i-1),2)+(K*pow(sin(1.0*u[i][j]),2)+1.0*pow(cos(1.0*u[i][j]),2))*pow(1.0*(i-1),4)*0.25*pow(v8+v9,2)+(K-1.0)*pow(1.0*(i-1),3)*0.5*(v8+v9)*sin(2.0*u[i][j]))
					-2.0*(K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1)*q[j],2)*(v1-v2)+4.0*(K*pow(sin(1.0*u[i][j]),2)+1.0*pow(cos(1.0*u[i][j]),2))*pow(1.0*(i-1)*q[j],2)+((K*sin(2.0*u[i][j])-sin(2.0*u[i][j]))*pow(1.0*(i-1),3)*(v8+v9)*q[j]+2.0*(K-1.0)*q[j]*pow(1.0*(i-1),2)*cos(2.0*u[i][j]))*0.5*(v6+v7);*/
		
}
