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

int main()
{
    int i,j,k,n,ncycle=2,ii,ii2,ii3,ii4,sos;
	double **u,x[300],y[300],*q,*dq,H,h,x1,x2,x3,fx1,fx2,fx3,xacc,R0,*vol,VO,s,rss;
	void mgfas(double **u, int n,int ncycle, double m, int a,double *q);
	void **nrfunc(double **f,double **u, int n);
	void rtsafe(double H, double h,double *L,int n,double **u,double *q,double *dq,int j,double R0,double t,double *su,double s,double rss);
	void rtsafe1( double x1, double x2, double H, double h,double *L,int n,double **u,double *q,double *dq, int j,double R0,double t,double*su,double s, double rss);
	/*void fanddf (double *f,double *df,double **u, double *q, double H, double h,double *L,int n,int j);*/
	void vofunc (double *q, double *vol,int n);
	double Stat (double **u,int n);
	/*void fanddf (double *f,double *df,double **u, double *q, double H, double h,double *L,int n,int j,double t,double s,int sos,double rss);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,double **u, double *q, double H, double h,double *L,int n,int j,double t,double s,double rss) ;
	*/double result,result2,**f,zero,rev2,rev3,v1,v2,v4,La,t,sum,result22,sum2,fff,dfff;
	double m; int a,b; 
	FILE *fp,*fp2,*fp3;
	n=65;
	
	La=/*-10843.421805224*/-0.1;xacc=0.00001;
	m=-1.0; a=32; b=50;R0=/*3.1942345*/10.0;rss=9.8784E14;printf("rss=%f\n",rss);
	h=1.0/(n-1); H=1.0;s=3500000000; 
	t=0.1;sum=0.0;

	u=dmatrix2(1,n,1,n);
	q=dvector(1,n);
	dq=dvector(1,n);
	
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
			for(k=1;k<=n;k++)
			    u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=a;i++)
		u[i][1]=Down1-0.5*m*Pi;
	
	for (i=1+a;i<=n;i++)
		u[i][n]=Down2;

	for (j=1;j<=n;j++)
		{
			u[n][j]=Right;
			u[1][j]=Left-0.5*m*Pi;
	    }
	

		mgfas(u,n,ncycle,m,a,q);

          j=33;
	                    v1=u[n][j+1]-u[n][j];
						v2=u[n][j]-u[n][j-1];
						v4=u[n][j]-u[n-1][j];

						if (v1>Pi/2+S)
						{v1=v1-Pi;}
						if (v1<-1.0*Pi/2-S)
							v1=v1+Pi;
						if (v2>Pi/2+S)
							v2=v2-Pi;
						if (v2<-1.0*Pi/2-S)
							v2=v2+Pi;
						if (v4>Pi/2+S)
							v4=v4-Pi;
						if (v4<-1.0*Pi/2-S)
							v4=v4+Pi;
		

La=/*(-pow(sin(u[n][j]),2)-v4*sin(u[n][j])/h+pow(sin(u[n][j]),2)*(v1+v2)/(h*H)-pow(v4*q[j]/h,2)-pow(q[j]*(v1+v2)/(2.0*h*H),2)*/(-A-B*(t*(q[j]-q[j-1])+(1.0-t)*(q[j+1]-q[j]))-C*pow(q[j],2)-D*pow((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1],2)-E*q[j]*((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1])-F)/pow(q[j],2);
printf("j=%d\tLa=%f\n",j,La);           

/*sum2=0.0;
for (j=2;j<=n-1;j++)
{
	                    v1=u[n][j+1]-u[n][j];
						v2=u[n][j]-u[n][j-1];
						v4=u[n][j]-u[n-1][j];

						if (v1>Pi/2+S)
						{v1=v1-Pi;}
						if (v1<-1.0*Pi/2-S)
							v1=v1+Pi;
						if (v2>Pi/2+S)
							v2=v2-Pi;
						if (v2<-1.0*Pi/2-S)
							v2=v2+Pi;
						if (v4>Pi/2+S)
							v4=v4-Pi;
						if (v4<-1.0*Pi/2-S)
							v4=v4+Pi;
						result22=abs(A+B*((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1])+C*pow(q[j],2)+D*pow((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1],2)+E*q[j]*((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1])+La*pow(q[j],2)+F);
						sum2=sum2+result22;
}

printf("La=%f\tsum2=%f\n",La,sum2);*/

/*fff=1.0;dfff=1.0;j=32;q[j]=11.0; t=0.0;s=0.0;
             fanddf (&fff,&dfff,u,q,H,h,&La,n,j,t,s);        
			 printf ("q[%d]=%f\n,f=%f\n,df=%f\n",j,q[j],fff,dfff);*/

			sum=0.0;t=0.0;
			for (j=32;j>=2;j--)
			{
                x1=0.0*R0;x2=1.1*R0;
				/*mnbrak(&x1, &x2, &x3, &fx1, &fx2, &fx3,u, q,  H,  h,&La,n,j, t, s);
				printf("x1=%f\nx2=%f\nx3=%f\nfx1=%f\nfx2=%f\nfx3=%f\n",x1,x2,x3,fx1,fx2,fx3);*/
				rtsafe1(x1, x2,H,h, &La,n,u,q,dq,j,R0,t,&sum, s,rss); 
				printf("Number 1 q[%d]=%f\n",j,q[j]);
			}

				/*for (j=32;j>=2;j--)
			{
               
				rtsafe( H, h,&La,n,u,q,dq,j,R0,t,&sum,s);
				printf("Number 2 q[%d]=%f\n",j,q[j]);
			}*/
			t=1.0;
			for (j=34;j<=n-1;j++)
			{
                x1=0.0*R0;x2=1.1*R0;
				rtsafe1(x1, x2,H,h, &La,n,u,q,dq,j,R0,t,&sum, s,rss); 
				printf("Number 1 q[%d]=%f\n",j,q[j]);
			}

			/*for (j=1;j<=n;j++)
				printf("q[%d]=%f\n",j,q[j]);*/

/*t=1.0;
		for (j=34;j<=n-1;j++)
	{ 
		        x1=9.5;x2=10.5;
				rtsafe1(x1, x2,H,h, &La,n,u,q,dq,j,R0, t,&sum, s); 
				printf("q[%d]=%f\n",j,q[j]);
}*/
printf("La=%f\tsum=%f\n",La,sum);

/*vofunc (q,vol,n);
	VO=trapzd(vol, n);
	printf("volume=%f\n",VO);
	for (j=1;j<=n;j++)
		q[j]=q[j]*sqrt(pow(R0,2)/VO);*/

/*for (ii3=1;ii3<=10;ii3++)
	{
		
		for (j=1;j<=b;j++)
		q[j]=q[33];
	for (j=b+1;j<=n;j++)
		q[j]=q[j-1];
	mgfas(u,n,ncycle,m,a,q);

                        j=33;
	                  v1=u[n][j+1]-u[n][j];
						v2=u[n][j]-u[n][j-1];
						v4=u[n][j]-u[n-1][j];

						if (v1>Pi/2+S)
						{v1=v1-Pi;}
						if (v1<-1.0*Pi/2-S)
							v1=v1+Pi;
						if (v2>Pi/2+S)
							v2=v2-Pi;
						if (v2<-1.0*Pi/2-S)
							v2=v2+Pi;
						if (v4>Pi/2+S)
							v4=v4-Pi;
						if (v4<-1.0*Pi/2-S)
							v4=v4+Pi;
		

La=(-A-B*(t*(q[j]-q[j-1])+(1.0-t)*(q[j+1]-q[j]))-C*pow(q[j],2)-D*pow((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1],2)-E*q[j]*((1.0-t)*q[j+1]+(2.0*t-1.0)*q[j]-t*q[j-1])-F)/pow(q[j],2);
printf("La=%f\n",La);

       sum=0.0;t=0.0;
			for (j=32;j>=2;j--)
			{
				rtsafe( H, h,&La,n,u,q,dq,j,R0,t,&sum,s);
                
			}
       t=1.0;
		for (j=34;j<=n-1;j++)
	{ 
		rtsafe( H, h,&La,n,u,q,dq,j,R0,t,&sum,s);

}

printf("La=%f\tsum=%f\n",La,sum);


vofunc (q,vol,n);
	VO=trapzd(vol, n);
	printf("volume=%f\n",VO);
	for (j=1;j<=n;j++)
		q[j]=q[j]*sqrt(pow(R0,2)/VO);

	}*/



/*for (ii2=1;ii2<=10;ii2++)
{
	                    j=33;
	                    v1=u[n][j+1]-u[n][j];
						v2=u[n][j]-u[n][j-1];
						v4=u[n][j]-u[n-1][j];

						if (v1>Pi/2+S)
						{v1=v1-Pi;}
						if (v1<-1.0*Pi/2-S)
							v1=v1+Pi;
						if (v2>Pi/2+S)
							v2=v2-Pi;
						if (v2<-1.0*Pi/2-S)
							v2=v2+Pi;
						if (v4>Pi/2+S)
							v4=v4-Pi;
						if (v4<-1.0*Pi/2-S)
							v4=v4+Pi;
		

La=(-pow(sin(u[n][j]),2)+pow(sin(u[n][j]),2)*(v1+v2)/(h*H)-pow(R0*(v1+v2)/(2.0*h*H),2))/pow(q[j],2);
printf("j=%d\tLa=%f\n",j,La);
			sum=0.0;			
		for (j=2;j<=n-1;j++)
	{ 
		rtsafe( H, h,&La,n,u,q,dq,j,R0,t,&sum);

}
printf("La=%f\tsum=%f\n",La,sum);

}*/

  printf("s=%f\n",s);                    

	 fp = fopen ("C:\\C\\NLPDE.txt","w");
	 fp2 = fopen ("C:\\C\\VectorField.txt","w");
	 fp3 = fopen ("C:\\C\\Shape.txt","w");
    if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
        return 1;
    }
	if (fp2 == NULL) {
        printf ("File 2 not created okay, errno = %d\n", errno);
        return 1;
    }
	if (fp3 == NULL) {
        printf ("File 3 not created okay, errno = %d\n", errno);
        return 1;
    }
		fprintf(fp,"{");
		for(i=1;i<=n-1;i++)
		{
			fprintf(fp,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp,"{%f,\t%f,\t%f},\t",x[i],y[k],u[i][k]);
			fprintf(fp,"{%f,\t%f,\t%f}\t",x[i],y[n],u[i][n]);
			fprintf(fp,"},");
		}
			fprintf(fp,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp,"{%f,\t%f,\t%f},\t",x[n],y[k],u[n][k]);
			fprintf(fp,"{%f,\t%f,\t%f}\t",x[n],y[n],u[n][n]);
			fprintf(fp,"}");
		fprintf(fp,"}");
		fclose (fp);
		fprintf(fp2,"{");
		for(i=1;i<=n-1;i++)
		{
			fprintf(fp2,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(u[i][k]),cos(u[i][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(u[i][n]),cos(u[i][n]));
			fprintf(fp2,"},");
		}
			fprintf(fp2,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(u[n][k]),cos(u[n][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(u[n][n]),cos(u[n][n]));
			fprintf(fp2,"}");
		fprintf(fp2,"}");
		fclose (fp2);

		fprintf(fp3,"{");
		for(j=1;j<=n-1;j++)
		{
			/*fprintf(fp3,"{");*/
			fprintf(fp3,"{%f,\t%f},\t",y[j],q[j]);
			/*fprintf(fp3,"},");*/
		}
		fprintf(fp3,"{%f,\t%f}\t",y[n],q[n]);
		fprintf(fp3,"}");
		fclose (fp3);
    printf ("File created okay\n");

	/*f=dmatrix2(1,n,1,n);
	nrfunc(f,u,n);
	zero=Stat(u,n);
	result=integrat(f, n);
	
	printf("%d\t%f\t%f\tf[n][n/2]=%f\nzero=%f\n",a,m,result,f[n][n/2],zero);*/

	/*m=-1.0*m;
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
			for(k=1;k<=a;k++)
			    u[i][k]=Down+0.5*m*Pi;
			for (k=1+a;k<=n;k++)
				u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=n;i++)
		u[i][1]=Down+0.5*m*Pi;
	
	for (i=1;i<=n;i++)
		u[i][n]=Up-0.5*m*Pi;
	for (j=2;j<=n/2;j++)
		u[n][j]=Right; 
	for (j=1+n/2;j<=n-1;j++)
		u[n][j]=Right; 
	mgfas(u,n,ncycle,m,a);
	nrfunc(f,u,n);
	zero=Stat(u,n);
	result2=integrat(f, n);
	if (result2<result)
	{
		result=result2;
		 fp = fopen ("C:\\C\\NLPDE.txt","w");
	 fp2 = fopen ("C:\\C\\VectorField.txt","w");
    if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
        return 1;
    }
	if (fp2 == NULL) {
        printf ("File 2 not created okay, errno = %d\n", errno);
        return 1;
    }
		
		fprintf(fp,"{");
		for(i=1;i<=n-1;i++)
		{
			fprintf(fp,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp,"{%f,\t%f,\t%f},\t",x[i],y[k],u[i][k]);
			fprintf(fp,"{%f,\t%f,\t%f}\t",x[i],y[n],u[i][n]);
			fprintf(fp,"},");
		}
			fprintf(fp,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp,"{%f,\t%f,\t%f},\t",x[n],y[k],u[n][k]);
			fprintf(fp,"{%f,\t%f,\t%f}\t",x[n],y[n],u[n][n]);
			fprintf(fp,"}");
		fprintf(fp,"}");
		fclose (fp);
		fprintf(fp2,"{");
		for(i=1;i<=n-1;i++)
		{
			fprintf(fp2,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(u[i][k]),cos(u[i][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(u[i][n]),cos(u[i][n]));
			fprintf(fp2,"},");
		}
			fprintf(fp2,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(u[n][k]),cos(u[n][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(u[n][n]),cos(u[n][n]));
			fprintf(fp2,"}");
		fprintf(fp2,"}");
		
		fclose (fp2);
    printf ("File created okay\n");
	}
 
	printf("%d\t%f\t%f\tf[n][n/2]=%f\nzero=%f\n",a,m,result2,f[n][n/2],zero);

	for (a=3;a<=n-1;a++)
	{
		m=-1.0*m;
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
			for(k=1;k<=a;k++)
			    u[i][k]=Down+0.5*m*Pi;
			for (k=1+a;k<=n;k++)
				u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=n;i++)
		u[i][1]=Down+0.5*m*Pi;
	
	for (i=1;i<=n;i++)
		u[i][n]=Up-0.5*m*Pi;
	for (j=2;j<=n/2;j++)
		u[n][j]=Right; 
	for (j=1+n/2;j<=n-1;j++)
		u[n][j]=Right; 
	mgfas(u,n,ncycle,m,a);
	nrfunc(f,u,n);
	zero=Stat(u,n);
	result2=integrat(f, n);
	if (result2<result)
	{
		result=result2;
		 fp = fopen ("C:\\C\\NLPDE.txt","w");
	 fp2 = fopen ("C:\\C\\VectorField.txt","w");
    if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
        return 1;
    }
	if (fp2 == NULL) {
        printf ("File 2 not created okay, errno = %d\n", errno);
        return 1;
    }
		
		fprintf(fp,"{");
		for(i=1;i<=n-1;i++)
		{
			fprintf(fp,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp,"{%f,\t%f,\t%f},\t",x[i],y[k],u[i][k]);
			fprintf(fp,"{%f,\t%f,\t%f}\t",x[i],y[n],u[i][n]);
			fprintf(fp,"},");
		}
			fprintf(fp,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp,"{%f,\t%f,\t%f},\t",x[n],y[k],u[n][k]);
			fprintf(fp,"{%f,\t%f,\t%f}\t",x[n],y[n],u[n][n]);
			fprintf(fp,"}");
		fprintf(fp,"}");
		fclose (fp);
		fprintf(fp2,"{");
		for(i=1;i<=n-1;i++)
		{
			fprintf(fp2,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(u[i][k]),cos(u[i][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(u[i][n]),cos(u[i][n]));
			fprintf(fp2,"},");
		}
			fprintf(fp2,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(u[n][k]),cos(u[n][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(u[n][n]),cos(u[n][n]));
			fprintf(fp2,"}");
		fprintf(fp2,"}");
		
		fclose (fp2);
    printf ("File created okay\n");
	}

	printf("%d\t%f\t%f\tf[n][n/2]=%f\nzero=%f\n",a,m,result2,f[n][n/2],zero);
		m=-1.0*m;
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
			for(k=1;k<=a;k++)
			    u[i][k]=Down+0.5*m*Pi;
			for (k=1+a;k<=n;k++)
				u[i][k]=Up-0.5*m*Pi;
	    }
	for (i=1;i<=n;i++)
		u[i][1]=Down+0.5*m*Pi;
	
	for (i=1;i<=n;i++)
		u[i][n]=Up-0.5*m*Pi;
	for (j=2;j<=n/2;j++)
		u[n][j]=Right; 
	for (j=1+n/2;j<=n-1;j++)
		u[n][j]=Right; 

		mgfas(u,n,ncycle,m,a);

	nrfunc(f,u,n);
	zero=Stat(u,n);
	result2=integrat(f, n);
	if (result2<result)
	{
		result=result2;
		 fp = fopen ("C:\\C\\NLPDE.txt","w");
	 fp2 = fopen ("C:\\C\\VectorField.txt","w");
    if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
        return 1;
    }
	if (fp2 == NULL) {
        printf ("File 2 not created okay, errno = %d\n", errno);
        return 1;
    }
		
		fprintf(fp,"{");
		for(i=1;i<=n-1;i++)
		{
			fprintf(fp,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp,"{%f,\t%f,\t%f},\t",x[i],y[k],u[i][k]);
			fprintf(fp,"{%f,\t%f,\t%f}\t",x[i],y[n],u[i][n]);
			fprintf(fp,"},");
		}
			fprintf(fp,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp,"{%f,\t%f,\t%f},\t",x[n],y[k],u[n][k]);
			fprintf(fp,"{%f,\t%f,\t%f}\t",x[n],y[n],u[n][n]);
			fprintf(fp,"}");
		fprintf(fp,"}");
		fclose (fp);
		fprintf(fp2,"{");
		for(i=1;i<=n-1;i++)
		{
			fprintf(fp2,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(u[i][k]),cos(u[i][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(u[i][n]),cos(u[i][n]));
			fprintf(fp2,"},");
		}
			fprintf(fp2,"{");
			for(k=1;k<=n-1;k++)
				fprintf(fp2,"{%f,\t%f},\t",sin(u[n][k]),cos(u[n][k]));
			fprintf(fp2,"{%f,\t%f}\t",sin(u[n][n]),cos(u[n][n]));
			fprintf(fp2,"}");
		fprintf(fp2,"}");
		
		fclose (fp2);
    printf ("File created okay\n");
	}
	
	printf("%d\t%f\t%f\tf[n][n/2]=%f\nzero=%f\n",a,m,result2,f[n][n/2],zero);
	}*/
	/*free_dmatrix2(f,1,n,1,n);*/
	free_dmatrix2(u,1,n,1,n);
	free_dvector(q,1,n);
	free_dvector(dq,1,n);
	free_dvector(vol,1,n);
    return 0;
}