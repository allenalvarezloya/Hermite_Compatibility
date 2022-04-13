#ifndef __HERMITE_H_INCLUDED__

#define __HERMITE_H_INCLUDED__
#include <math.h>
#include <iostream>
#include "Darrays.h"
#include "extern.h"
extern "C" void tpolymul2_(double *,double *,double *, int *);
extern "C" void dgetrf_(int *, int *, double *, int *, int *, int *);
extern "C" void dgemv_(char *,int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
extern "C" void dgetrs_(char *,int *, int *, double *, int *, int *, double *, int *, int *);
extern "C" void dgesvd_(char *,char *,int *,int *,double *,int *,double *,double *,int *,double *,int *,double *,int *,int *);
extern "C" void eq_(double *,int *,int *,int *,double *,double *,double *);
class Hermite {
public:
    int m, nr, ns;
    Darray2 Hmat_u, Hmat_v;
    void binomial(int m,Darray2 &coeffs);
    void HermiteMap(int m,double xl,double xr,double xc,int icase,Darray2 &tmat);
    void compute_Hmaps();
    void interpPD(Darray2 &u,Darray2 &v);
    void interpDP(Darray2 &ud,Darray2 &vd);
	Darray2 u,u_interp,ud,ud_interp,v,v_interp,vd,vd_interp;
    void outPutVec(Darray1 &M,int i_start,int i_end,char* fName);
    void outPutMat(Darray2 &M,int i_start,int i_end,int j_start,int j_end,char* fName);
    void recursion(Darray2 &u_int, Darray2 &v_int,Darray2 &u, Darray2 &v,int i,double h,double C,double dt);
    void bcMat(double rA,double hr,double *MM,int m_mat);
    void boundaryConditions(int m, double hr,Darray2 &u, Darray2 &v);
    void interface(double hu,double hv,double C,double D,double *MM,int m_mat);
    void interfaceConditions(int m,double hu,double hv,double C,double D,Darray2 &u1,Darray2 &u2, Darray2 &v1, Darray2 &v2,int i1,int i2,int iout,int side);
    void oversample(Darray2 &u,int nr,int m,char* fName);
};

#endif