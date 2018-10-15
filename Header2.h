#ifndef HEADER2_H
#define HEADER2_H
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix2(long nrl,long nrh,long ncl,long nch);
void free_dmatrix2(double **m, long nrl, long nrh, long ncl, long nch);
void fill02(double **u, int n, double m, int a);
void free_vector(float *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix3(double **m, long nrl, long nrh, long ncl, long nch);
void free_ivector(int *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
#endif