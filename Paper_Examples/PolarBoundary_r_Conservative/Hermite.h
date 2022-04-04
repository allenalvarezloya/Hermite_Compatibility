#ifndef __HERMITE_H_INCLUDED__

#define __HERMITE_H_INCLUDED__
#include <math.h>
#include <iostream>
#include "Darrays.h"
#include "extern.h"
#include <boost/math/special_functions/bessel.hpp>
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
    void interpolate_pd(Darray4 &u); // Primal to Dual
	void interpolate_dp(Darray4 &ud); //  Dual to primal
	Darray4 u,u_interp,ud,ud_interp,v,v_interp,vd,vd_interp;
    void cartesian_recursion(Darray4 &u_int,Darray4 &v_int,Darray4 &u,Darray4 &v,int i,int j,double hr, double hs, double dt);
	void recursion_curvilinear(Darray4 &u_int,int i,int j,
    Darray4 &A10, Darray4 &A01,Darray4 & A11,Darray4 &A12,Darray4 & A22,Darray4 &u,
    double hr,double hs,double dt);
	void bcMat(double rA,double rB,double sA,double sB,double hr,double hs,
        double *a01,double* a10,double* a11,double* a12,double* a22,double *MM,int m_mat);
    void boundaryConditions(int m,double hr,double hs,Darray4 &u,Darray4 &A01,Darray4 & A10,Darray4 & A11,Darray4 &A12,Darray4 & A22);
    void BoundaryPeriodic(Darray4 &ud);
    void oversample(Darray4 &u,Darray2 &U,int nsample_r,int nsample_s);
    double HBessel(double nu,double x);
    void outPutMat(Darray2 &M,int i_start,int i_end,int j_start,int j_end,char* fName);
    void outPutVec(Darray1 &M,int i_start,int i_end,char* fName);
    void computeSingular(int m,Darray2 Mat_u,int ns);
    void curviLaplacian(double hx,double hy,int m,double dt,Darray2 &a01,Darray2 &a10,Darray2 &a11,Darray2 &a12,Darray2 &a22,Darray2 &u);
};

#endif