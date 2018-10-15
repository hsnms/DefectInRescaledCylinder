#ifndef HEADER3_H
#define HEADER3_H
double anorm2(double **a, int n);
void matadd(double **a, double **b, double **c, int n);
void matsub(double **a, double **b, double **c, int n);
void copy2(double **aout, double **ain, int n);
void relax2(double **u, double **rhs, int n,double m,int a,double *q);
void interp2(double **uf, double **uc, int nf, double m, int a); 
void interp3(double **uf, double **uc, int nf, int a);
void rstrct2(double **uc, double **uf, int nc);
void lop(double **out, double **u, int n, double m, int a,double *q);
void vecfunc(int n, int k, float *x, float *f, double m, int a);
#endif