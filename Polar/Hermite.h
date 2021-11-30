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
class Hermite {
public:
    int m, nr, ns;
    Darray2 Hmat_u, Hmat_v;
    void binomial(int m,Darray2 &coeffs);
    void HermiteMap(int m,double xl,double xr,double xc,int icase,Darray2 &tmat);
    void compute_Hmaps();
    void interpolate_pd(Darray4 &u, Darray4 &v); // Primal to Dual
	void interpolate_dp(Darray4 &ud, Darray4 &vd); //  Dual to primal
    void interpolate_dp_ext(Darray4 &ud, Darray4 &vd); //  Dual to primal
	Darray4 u,u_interp,ud,ud_interp,v,v_interp,vd,vd_interp;
    void cartesian_recursion(Darray4 &u_int,Darray4 &v_int,Darray4 &u,Darray4 &v,int i,int j,double hr, double hs, double dt);
	void recursion_curvilinear(Darray4 &u_int,Darray4 &v_int,int i,int j,
    Darray4 &A10, Darray4 &A01,Darray4 & A11,Darray4 &A12,Darray4 & A22,Darray4 &u,Darray4 &v,
    double hr,double hs,double dt);
	void M1(double rA,double rB,double sA,double sB,double hr,double hs,
        double *a01,double* a10,double* a11,double* a12,double* a22,double *MM);
    void M2(double rA,double rB,double sA,double sB,double hr,double hs,
        double *a01,double* a10,double* a11,double* a12,double* a22,double *MM);
    // void M3(double rA,double rB,double sA,double sB,double hr,double hs,
    //     double *a01,double* a10,double* a11,double* a12,double* a22,double *MM);
    void boundaryConditions(int m,double hr,double hs,Darray4 &u,Darray4 &v,Darray4 &A01,Darray4 & A10,Darray4 & A11,Darray4 &A12,Darray4 & A22);
    void BoundaryPeriodic(Darray4 &ud, Darray4 &vd);
    void oversample(Darray4 &u,double time,Darray1 &errors1,Darray1 &errors2,Darray1 &errorsinf,int COUNT);
    double HBessel(double nu,double x);
};

#endif