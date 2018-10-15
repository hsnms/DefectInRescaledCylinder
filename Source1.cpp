//Source1.cpp helps to find the specific equilibrium nematic defect structures in a cylinder given specific parameters including the Frank constant ratio and aspect ratio.

//We consider homeotropic boundary conditions and proper inner boundary conditions (after we apply a cut-off of the defect core), and solve the Euler-Lagrange equation. Because of the cylindrical symmetry, we only consider the rectangular diametrical plane.

//The finite difference method (particularly the over-relaxation method) is used. The rectangular region is rescaled to a square, the aspect ratio appears in the equation and boundary condition.

//Input: n: number of lattices horizontally and vertically; R0: aspect ratio; K:Frank constant ratio; b: the number of the lattice spaces the length of the defect core occupies; m:types of defect structures (radial:1.0, hyperbolic:-1.0).

//Output: 1,"Bridge.txt": different energies when the defect core is at different locations, error, all the input values; 2, "BridgeVectorField.txt": the vector field.

//How to compile: 1, make sure the in the "makefile", we have the relevant "Source1"; 2, M-X compile (option + X + compile); 3, Compile command: make -k output; 4, Save file, y; 5, M-X shell; 6, bash-3.2$ ./output









#include <stdio.h>
#include <errno.h>
#include <math.h>
#define Up 1.570796326794897
#define Down1 1.570796326794897
#define Down2 1.570796326794897
#define Right 1.570796326794897
#define Left 1.570796326794897
#define Pi 3.141592653589793
#define S 0.00000000000001
#include "Header1.h"
#include "Header2.h"
#include "Header5.h"
#include "Header6.h"
int main()
{
    int i,j,k,n,ncycle=5,ii,ii2,ii3,ii4,sos,kk,radius,min,n1;
	double **u,*x,*y,*z,*z2,*q,*dq,H,h,x1,x2,x3,fx1,fx2,fx3,xacc,R0,*vol,VO,s,rss,**v;
	void mgfas4(double **u, int n,int ncycle, double m, int a,int b,double *q, double *x, double *y, double R0, double K);
	void mgfas5(double **u, int n,int ncycle, double m,int b,double *q, double *x, double *y, double R0, double K);
	void **nrfunc2(double **f,double **u, int n,double *q, int a);
	void **nrfunc3(double **f,double **u, int n,double *q, int a,int b,double R0, double K);
	void **nrfunc4(double **f,double **u, int n,double *q,int b,double R0, double K); 
	double Stat2 (double **u,int n,double *q, int a,int b, double R0, double K);
	double Stat3 (double **u,int n,double *q,int b, double R0, double K);
	double result,/*result2,*/**f,zero,rev2,rev3,v1,v2,v4,K,t,sum,result22,sum2,fff,dfff;
	double m,point,point2; int a,b; 
	FILE *fp,*fp2,*fp3;

	printf ("Number of lattices horizontally or vertically, n:");
	scanf("%d", &n);
	//n=17 (or 33, 65, 129, 257);
        printf ("Aspect ratio, R0:");
	scanf("%lf", &R0);
	//R0=1;
        H=1.0;
        printf ("Frank constant ratio, K:");
	scanf("%lf", &K);
	//K=1;
 
	fp = fopen ("Bridge.txt","w");
	  if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
        return 1;
    }
	  printf ("How many of the lattice spaces does the length of the defect core occupy? b:");
	scanf("%d", &b);
	// b=1 (or 2, 4, 8);
	kk=1+2*b;
        printf ("Types of defect structures (radial:1.0, hyperbolic:-1.0), m:");
	scanf("%lf", &m);
	
	//m=-1.0;

	a=kk;//smallest radius of the ring defect
    u=dmatrix2(1,n,1,n);
	v=dmatrix2(1,n,1,n);
	q=dvector(1,n);
	x=dvector(1,n);
	y=dvector(1,n);
	z=dvector(1,n);
	z2=dvector(1,n);
	for (i=1;i<=n;i++)//initial and outer boundary conditions
	{
		x[i]=1.0*(i-1)*(1.0/(n-1));
	}
	for (i=1;i<=n;i++)
	{
		y[i]=1.0*(i-1)*(1.0/(n-1));
	} 
        for(i=1;i<=n;i++)
       	{
			    u[i][n]=Up-0.5*m*Pi;
	 }
	for(i=2;i<=n-1;i++)
		{
			for(k=2;k<=n-1;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=b;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+b;i<=n;i++)
		u[i][1]=Down2;

	for (j=1;j<=n;j++)
		{
			u[n][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }

	/***/
	if (m==-1.0)// hyperbolic type
	{

	  //point defect
	  //inner boundary conditions
	for (j=2;j<=1+b*(int)(R0/H);j++)
	{
	  //u[1+b][j]=atan(H*(y[j]-y[1])/(R0*(x[1+b]-x[1])))+0.5*Pi-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*atan(H*(y[j]-y[1])/(R0*(x[1+b]-x[1])))));
	  u[1+b][j]=atan(H*(y[j]-y[1])/(R0*(x[1+b]-x[1])))+0.5*Pi-((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*atan(H*(y[j]-y[1])/(R0*(x[1+b]-x[1])))));
		
	}
	for (i=1;i<=1+b;i++)
	{
	  //	u[i][1+b*(int)(R0/H)]=atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[1])))+0.5*Pi-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[1])))));
	  u[i][1+b*(int)(R0/H)]=atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[1])))+0.5*Pi-((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[1])))));	
	}
	
	for (i=n-b;i<=n-1;i++)
	{
		u[i][n-b*(int)(R0/H)]=Pi-atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])));
		
	}
	for (j=n-b*(int)(R0/H);j<=n-1;j++)
    {
		u[n-b][j]=Pi-atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])));
		
	}
	for (j=1;j<=n;j++)
		q[j]=R0;
	mgfas5(u,n,ncycle,m,b,q,x,y, R0,K);//solve PDE
		f=dmatrix2(1,n,1,n);
		if (R0>=1.0)
		{
		  nrfunc4(f,u,n,q,b,R0,K);//integrand
		}
		else
		{
			nrfunc2(f,u,n,q,kk);
		}
n1=1;
 z[1]=integrat(f, n1,n,n1,n,n)-integrat(f, n1,n1+b,n1,b*(int)(R0/H)+1,n)-integrat(f, n-b,n,n-b*(int)(R0/H),n,n);//total energy
 zero=Stat3(u,n,q,b,R0,K);//error
	result=z[1];
	for(i=1;i<=n;i++)
		for(j=1;j<=n;j++)
			v[i][j]=u[i][j];
	radius=1;
	
	printf("z[1]=%f\t\nzero=%f\n",z[1],zero);
	fprintf(fp,"z[1]=%f\t\n{0, %f},\nzero=%f\n",z[1],z[1],zero);


	//ring defect with smallest radius
	//initial and outer boundary conditions
		for(i=1;i<=n;i++)
		{
			    u[i][n]=Up-0.5*m*Pi;
	    }
	for(i=2;i<=n-1;i++)
		{
			for(k=2;k<=n-1;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=a;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=n;i++)
		u[i][1]=Down2;

	for (j=1;j<=n;j++)
		{
			u[n][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }

	//inner boundary conditions
	for (j=2;j<=1+b*(int)(R0/H);j++)
	{
	  //	u[a+b][j]=0.5*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))));
	  u[a+b][j]=0.5*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))));	  
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b*(int)(R0/H)]=0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))));
	  u[i][1+b*(int)(R0/H)]=0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))));
	}
//	u[a][1+b*(int)(R0/H)]=0.5*0.5*Pi+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
u[a][1+b*(int)(R0/H)]=0.5*0.5*Pi+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
      //	u[i][1+b*(int)(R0/H)]=-0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i]))))));
      u[i][1+b*(int)(R0/H)]=-0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i]))))));
	}
	for (j=2;j<=1+b*(int)(R0/H);j++)
    {
      //	u[a-b][j]=-0.5*atan(H*(y[j]-y[1])/(R0*(x[a]-x[a-b])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[a-b])))))); the original one is wrong
      u[a-b][j]=-0.5*atan(H*(y[j]-y[1])/(R0*(x[a]-x[a-b])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan(H*(y[j]-y[1])/(R0*(x[a]-x[a-b]))))));
	}
	for (i=n-b;i<=n-1;i++)
	{
		u[i][n-b*(int)(R0/H)]=Pi-atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])));
	}
	for (j=n-b*(int)(R0/H);j<=n-1;j++)
    {
		u[n-b][j]=Pi-atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])));
	}
	for (j=1;j<=n;j++)
		q[j]=R0;
        mgfas4(u,n,ncycle,m,a,b,q,x,y, R0,K); //solve PDE
		f=dmatrix2(1,n,1,n);
		if (R0>=1.0)
		{
		  nrfunc3(f,u,n,q,kk,b,R0,K);//integrand
		}
		else
		{
			nrfunc2(f,u,n,q,kk);
		}
n1=1;
 z[kk]=integrat(f, n1,n,n1,n,n)-integrat(f, a-b,a+b,n1,b*(int)(R0/H)+1,n)-integrat(f, n-b,n,n-b*(int)(R0/H),n,n);//total energy 
 zero=Stat2(u,n,q,a,b,R0,K);//error
  
 //find energy minimum
        if (z[kk]<result)
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=n;i++)
		for(j=1;j<=n;j++)
			v[i][j]=u[i][j];
	}

	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);


	
	//ring defect with different radii
	//initial and outer boundary conditions
	for(kk=1+2*b;kk<=n-b-1;kk++)
	{
	 a=kk;
	

	u=dmatrix2(1,n,1,n);
	q=dvector(1,n);
	
	
	for (i=1;i<=n;i++)
	{
		x[i]=1.0*(i-1)*(1.0/(n-1));
	}
	for (i=1;i<=n;i++)
	{
		y[i]=1.0*(i-1)*(1.0/(n-1));
	}




for(i=1;i<=n;i++)
		{
			    u[i][n]=Up-0.5*m*Pi;
	    }
	for(i=2;i<=n-1;i++)
		{
			for(k=2;k<=n-1;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=a;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=n;i++)
		u[i][1]=Down2;

	for (j=1;j<=n;j++)
		{
			u[n][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }

	//inner boundary conditions
	for (j=2;j<=1+b*(int)(R0/H);j++)
	{
	  //	u[a+b][j]=0.5*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))));
	  u[a+b][j]=0.5*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))));
	  
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b*(int)(R0/H)]=0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))));
	  u[i][1+b*(int)(R0/H)]=0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))));
       	}
//	u[a][1+b*(int)(R0/H)]=0.5*0.5*Pi+0.5*Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
	u[a][1+b*(int)(R0/H)]=0.5*0.5*Pi+0.5*Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
      //	u[i][1+b*(int)(R0/H)]=-0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i]))))));
      u[i][1+b*(int)(R0/H)]=-0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i]))))));
      
	}
	for (j=2;j<=1+b*(int)(R0/H);j++)
    {
      //	u[a-b][j]=-0.5*atan(H*(y[j]-y[1])/(R0*(x[a]-x[a-b])))+Pi-((1.0/(24.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[a-b]))))));
      u[a-b][j]=-0.5*atan(H*(y[j]-y[1])/(R0*(x[a]-x[a-b])))+Pi-((5.0/36.0)*((K-1.0)/(K+1.0))*sin(3.0*(Pi-atan(H*(y[j]-y[1])/(R0*(x[a]-x[a-b]))))));
	}
	for (i=n-b;i<=n-1;i++)
	{
		u[i][n-b*(int)(R0/H)]=Pi-atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])));
	}
	for (j=n-b*(int)(R0/H);j<=n-1;j++)
    {
		u[n-b][j]=Pi-atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])));
	}

	for (j=1;j<=n;j++)
		q[j]=R0;
        mgfas4(u,n,ncycle,m,a,b,q,x,y, R0,K);//solve PDE
		f=dmatrix2(1,n,1,n);
		if (R0>=1.0)
		{
		  nrfunc3(f,u,n,q,kk,b,R0,K);//integrand
		}
		else
		{
			nrfunc2(f,u,n,q,kk);
		}
	n1=1;
	z[kk]=integrat(f, n1,n,n1,n,n)-integrat(f, a-b,a+b,n1,b*(int)(R0/H)+1,n)-integrat(f, n-b,n,n-b*(int)(R0/H),n,n);//total energy
	zero=Stat2(u,n,q,a,b,R0,K); //error
	if (z[kk]<result)//find energy minimum
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=n;i++)
		for(j=1;j<=n;j++)
			v[i][j]=u[i][j];
	}
	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	
	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);
	}
}



	/***/
	//radial type
	if (m==1.0)
	{
	  //point defect
	  //inner boundary conditions
		for (j=2;j<=1+b*(int)(R0/H);j++)
	{
		u[1+b][j]=-atan(H*(y[j]-y[1])/(R0*(x[1+b]-x[1])))+0.5*Pi;
		
	}
	for (i=1;i<=1+b;i++)
	{
		u[i][1+b*(int)(R0/H)]=-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[1])))+0.5*Pi;
		
	}
	
	for (i=n-b;i<=n-1;i++)
	{
	  //	u[i][n-b*(int)(R0/H)]=atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))-Pi)));
	  u[i][n-b*(int)(R0/H)]=atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))-Pi)));
		
	}
	for (j=n-b*(int)(R0/H);j<=n-1;j++)
    {
      //	u[n-b][j]=atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))-Pi)));
      	u[n-b][j]=atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))-Pi)));
		
	}
	for (j=1;j<=n;j++)
		q[j]=R0;
	mgfas5(u,n,ncycle,m,b,q,x,y, R0,K);//solve PDE
		f=dmatrix2(1,n,1,n);
		if (R0>=1.0)
		{
		  nrfunc4(f,u,n,q,b,R0,K);//integrand
		}
		else
		{
			nrfunc2(f,u,n,q,kk);
		}
n1=1;
 z[1]=integrat(f, n1,n,n1,n,n)-integrat(f, n1,n1+b,n1,b*(int)(R0/H)+1,n)-integrat(f, n-b,n,n-b*(int)(R0/H),n,n);//total energy
 zero=Stat3(u,n,q,b,R0,K);//error
	result=z[1];
	for(i=1;i<=n;i++)
		for(j=1;j<=n;j++)
			v[i][j]=u[i][j];
	radius=1;
	
	printf("z[1]=%f\t\nzero=%f\n",z[1],zero);
	fprintf(fp,"z[1]=%f\t\n{0, %f},\nzero=%f\n",z[1],z[1],zero);


	//ring defect with smallest radius
	//initial and outer boundary conditions
		for(i=1;i<=n;i++)
		{
			    u[i][n]=Up-0.5*m*Pi;
	    }
	for(i=2;i<=n-1;i++)
		{
			for(k=2;k<=n-1;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=a;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=n;i++)
		u[i][1]=Down2;

	for (j=1;j<=n;j++)
		{
			u[n][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }
	//inner boundary conditions
for (j=2;j<=1+b*(int)(R0/H);j++)
	{
	  //	u[a+b][j]=-0.5*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))));
	  u[a+b][j]=-0.5*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))));
	  
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b*(int)(R0/H)]=-0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))));

	  	u[i][1+b*(int)(R0/H)]=-0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))));
	}
	u[a][1+b*(int)(R0/H)]=-0.5*0.5*Pi+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
      //	u[i][1+b*(int)(R0/H)]=0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin((Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i]))))));

      u[i][1+b*(int)(R0/H)]=0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin((Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i]))))));
	}
	for (j=2;j<=1+b*(int)(R0/H);j++)
    {
      //	u[a-b][j]=0.5*atan(H*(y[j]-y[1])/(R0*(x[a]-x[a-b])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[a-b]))))));

      	u[a-b][j]=0.5*atan(H*(y[j]-y[1])/(R0*(x[a]-x[a-b])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[a-b]))))));
	}
	for (i=n-b;i<=n-1;i++)
	{
	  //	u[i][n-b*(int)(R0/H)]=atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))-Pi)));
	  	  u[i][n-b*(int)(R0/H)]=atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))-Pi)));
		
	}
	for (j=n-b*(int)(R0/H);j<=n-1;j++)
    {
      //	u[n-b][j]=atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))-Pi)));
        	u[n-b][j]=atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))-Pi)));
		
	}
	for (j=1;j<=n;j++)
		q[j]=R0;
        mgfas4(u,n,ncycle,m,a,b,q,x,y, R0,K);//solve PDE
		f=dmatrix2(1,n,1,n);
		if (R0>=1.0)
		{
		  nrfunc3(f,u,n,q,kk,b,R0,K);//integrand
		}
		else
		{
			nrfunc2(f,u,n,q,kk);
		}
n1=1;
 z[kk]=integrat(f, n1,n,n1,n,n)-integrat(f, a-b,a+b,n1,b*(int)(R0/H)+1,n)-integrat(f, n-b,n,n-b*(int)(R0/H),n,n);//total energy
 zero=Stat2(u,n,q,a,b,R0,K);//error
  
 //find energy minimum
	if (z[kk]<result)
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=n;i++)
		for(j=1;j<=n;j++)
			v[i][j]=u[i][j];
	}

	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);


	//ring defect with different radii
	
	for(kk=1+2*b;kk<=n-b-1;kk++)
	{
	 a=kk;
	

	u=dmatrix2(1,n,1,n);
	q=dvector(1,n);
	
	//initial and outer boundary conditions
	for (i=1;i<=n;i++)
	{
		x[i]=1.0*(i-1)*(1.0/(n-1));
	}
	for (i=1;i<=n;i++)
	{
		y[i]=1.0*(i-1)*(1.0/(n-1));
	}




for(i=1;i<=n;i++)
		{
			    u[i][n]=Up-0.5*m*Pi;
	    }
	for(i=2;i<=n-1;i++)
		{
			for(k=2;k<=n-1;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=a;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=n;i++)
		u[i][1]=Down2;

	for (j=1;j<=n;j++)
		{
			u[n][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }

	//inner boundary conditions
	for (j=2;j<=1+b*(int)(R0/H);j++)
	{
	  //	u[a+b][j]=-0.5*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))));
  u[a+b][j]=-0.5*atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan(H*(y[j]-y[1])/(R0*(x[a+b]-x[a])))));
	  
	}
	for (i=a+1;i<=a+b;i++)
	{
	  //	u[i][1+b*(int)(R0/H)]=-0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))));

	  	u[i][1+b*(int)(R0/H)]=-0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[i]-x[a])))));
	}
	//	u[a][1+b*(int)(R0/H)]=-0.5*0.5*Pi+0.5*Pi+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(0.5*Pi));

		u[a][1+b*(int)(R0/H)]=-0.5*0.5*Pi+0.5*Pi+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(0.5*Pi));
	for (i=a-b;i<=a-1;i++)
    {
      //	u[i][1+b*(int)(R0/H)]=0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin((Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i]))))));

          u[i][1+b*(int)(R0/H)]=0.5*atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin((Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[i]))))));
	}
	for (j=2;j<=1+b*(int)(R0/H);j++)
    {
      //	u[a-b][j]=0.5*atan(H*(y[j]-y[1])/(R0*(x[a]-x[a-b])))+((5.0/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[a-b]))))));

       	u[a-b][j]=0.5*atan(H*(y[j]-y[1])/(R0*(x[a]-x[a-b])))+((3.0/4.0)*((K-1.0)/(K+1.0))*sin(1.0*(Pi-atan(H*(y[1+b*(int)(R0/H)]-y[1])/(R0*(x[a]-x[a-b]))))));
	}
	for (i=n-b;i<=n-1;i++)
	{
	  //	u[i][n-b*(int)(R0/H)]=atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))-Pi)));
	  	  u[i][n-b*(int)(R0/H)]=atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[n-b*(int)(R0/H)]-y[n])/(R0*(x[i]-x[n])))-Pi)));
		
	}
	for (j=n-b*(int)(R0/H);j<=n-1;j++)
    {
      //	u[n-b][j]=atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))-((0.25/(8.0*Pi*Pi))*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))-Pi)));
		u[n-b][j]=atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))+((3.0/16.0)*((K-1.0)/(K+1.0))*sin(4.0*(atan(H*(y[j]-y[n])/(R0*(x[n-b]-x[n])))-Pi)));	
	}
	for (j=1;j<=n;j++)
		q[j]=R0;
        mgfas4(u,n,ncycle,m,a,b,q,x,y, R0,K);//solve PDE
		f=dmatrix2(1,n,1,n);
		if (R0>=1.0)
		{
		  nrfunc3(f,u,n,q,kk,b,R0,K);//integrand
		}
		else
		{
			nrfunc2(f,u,n,q,kk);
		}
	n1=1;
	z[kk]=integrat(f, n1,n,n1,n,n)-integrat(f, a-b,a+b,n1,b*(int)(R0/H)+1,n)-integrat(f, n-b,n,n-b*(int)(R0/H),n,n);//total energy
	zero=Stat2(u,n,q,a,b,R0,K);//error
	if (z[kk]<result)//find energy minimum
	{
	result=z[kk];
	radius=kk;
	for(i=1;i<=n;i++)
		for(j=1;j<=n;j++)
			v[i][j]=u[i][j];
	}
	printf("z[%d]=%f\t\nzero=%f\n",kk,z[kk],zero);
	
	fprintf(fp,"z[%d]=%f\t\n{%f, %f},\nzero=%f\n",kk,z[kk],1.0*(kk-1)/(n-1),z[kk],zero);
	}
	}


			
	printf("Total energy=%f\t\nRing radius=%d\n",result,radius);
	fprintf(fp,"Number of lattices horizontally or vertically n=%d\nNumber of the lattice spaces the length of the core occupies b=%d\nDefect type m=%f\nAspect ratio R0=%f\nFrank constant ratio K=%f\nTotal energy=%f\nRing radius=%d\n",n,b,m,R0,K,result,radius);
	//fprintf(fp,"Number of lattices horizontally or vertically n=%d\nNumber of the lattice spaces the length of the core occupies b=%d\nDefect type m=%f\nAspect ratio R0=%f\nFrank constant ratio K=%f\nTotal energy=%f\nRing radius=%d\n",n,b,m,R0,K,result,radius);

	fclose (fp);
                   

	
	/* fp2 = fopen ("BridgeVectorField.txt","w");
  
	if (fp2 == NULL) {
        printf ("File 2 not created okay, errno = %d\n", errno);
        return 1;
    }
       
		
		

		fprintf(fp2,"{");
		for(i=1;i<=n-1;i++)
		{
			fprintf(fp2,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(v[i][k]),cos(v[i][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(v[i][n]),cos(v[i][n]));
			fprintf(fp2,"},");
		}
			fprintf(fp2,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(v[n][k]),cos(v[n][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(v[n][n]),cos(v[n][n]));
			fprintf(fp2,"}");
		fprintf(fp2,"}");
		fclose (fp2);*/
	
    printf ("File created okay\n");

	
	free_dmatrix2(u,1,n,1,n);
	free_dmatrix2(v,1,n,1,n);
	free_dvector(q,1,n);
	free_dvector(x,1,n);
	free_dvector(y,1,n);
	free_dvector(z,1,n);
	free_dvector(z2,1,n);
    return 0;
}
