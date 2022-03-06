#include <iostream>
#include <cmath>
#include "SetUp.h"
#include "Hermite.h"
using namespace std;
extern "C" void dgetrs_(char *,int *, int *, double *, int *, int *, double *, int *, int *);
int main(){
	int m = 2;
	int COUNT = 1;
	int plotCOUNT = 1;
	int refine = 1;
	Darray1 errors1, errors2, errorsinf;
	errors1.define(1,refine);
	errors2.define(1,refine);
	errorsinf.define(1,refine);
	errors1.set_value(0.0);
	errors2.set_value(0.0);
	errorsinf.set_value(0.0); 
	char fName[100];
	int base = 5;
	// For refinement
    int nsample_r, nsample_s;
    nsample_r = 10;    // Oversample in r
    nsample_s = 10;    // Oversample in s
	for(int n = base; n <= base*pow(2,refine); n *= 2){
	    double k31, k33,gamma,alpha;
	    k31 = 6.3801618959239383506237;
	    k33 = 13.01520072169843441983;
	    gamma = k31/k33;
	    alpha = 1-gamma;
		int nr = n;
		// double hr_phys = gamma/nr;
		// double ns = floor(2*M_PI*gamma/hr_phys); // Set ns
		// cout << "Number of cells in s = " << ns << endl;
		int ns = 5*n;


		// Allocate memory for solution arrays (Primal)
		Darray4 u,v;
		u.define(0,m,0,m,0,nr,0,ns);
		v.define(0,m-1,0,m-1,0,nr,0,ns);
		u.set_value(0.0);
		v.set_value(0.0);
		Darray2 U_refined;
		U_refined.define(0,nsample_r*nr,0,nsample_s*ns);
    	U_refined.set_value(0.0);
		Darray2 X_refined;
		X_refined.define(0,nsample_r*nr,0,nsample_s*ns);
    	X_refined.set_value(0.0);
		Darray2 Y_refined;
		Y_refined.define(0,nsample_r*nr,0,nsample_s*ns);
    	Y_refined.set_value(0.0);
		// Allocate memory for solution arrays (Dual)
		Darray4 ud,vd;
		ud.define(0,m,0,m,1,nr,1,ns);
		vd.define(0,m-1,0,m-1,1,nr,1,ns);
		ud.set_value(0.0);
		vd.set_value(0.0);
		// Setup information
		SetUp SetUp;
		SetUp.metricOrder = 2*m+5;
		SetUp.m = m;
		SetUp.nr = nr;
		SetUp.ns = ns;
		SetUp.x_s = 0.0;
		SetUp.x_e = 1.0;
		SetUp.y_s = 0.0;
		SetUp.y_e = 1.0;
		SetUp.compute_grids();
		SetUp.transformation();
		SetUp.compute_initial_data(u,v);
		SetUp.Metric_Info();
		SetUp.metricCoeffs();


		// Hermite information
		Hermite Hermite;
		Hermite.m = m;
		Hermite.nr = nr;
		Hermite.ns = ns;
		Hermite.compute_Hmaps();
		Hermite.u_interp.define(0,2*m+1,0,2*m+1,0,nr,0,ns);
		Hermite.u_interp.set_value(0.0);
		Hermite.v_interp.define(0,2*m-1,0,2*m-1,0,nr,0,ns);
		Hermite.v_interp.set_value(0.0);

		Darray2 xgrid, ygrid;
		xgrid.define(0,nr,0,ns);
		ygrid.define(0,nr,0,ns);
		xgrid.set_value(0.0);
		ygrid.set_value(0.0);
		for(int j = 0; j<=ns; j++){
			for(int i = 0; i <=nr; i++){
				xgrid(i,j) = SetUp.X(0,0,i,j);
				ygrid(i,j) = SetUp.Y(0,0,i,j);
			}
		}
    	sprintf(fName, "xgrid%05d.ext",COUNT);
    	Hermite.outPutMat(xgrid,0,nr,0,ns,fName);
    	sprintf(fName, "ygrid%05d.ext",COUNT);
    	Hermite.outPutMat(ygrid,0,nr,0,ns,fName);
  //   	double heff = 1.0e10;
		// for(int j = 1; j<=ns; j++){
		// 	for(int i = 1; i <=nr; i++){
		// 		double htemp =  min(abs(xgrid(i,j) - xgrid(i-1,j)),abs(ygrid(i,j)-ygrid(i,j-1)));
		// 		heff = (htemp < heff) ? htemp : heff;
		// 	}
		// }
		double heff = SetUp.computeMinJacobian(SetUp.hr,SetUp.hs,nr,ns,SetUp.X,SetUp.Y);
		double t, dt, T;
		double CFL = 0.5;
	    T = 2*M_PI/k33;
	    t = 0.0;
	    while (t < T){
	    	if (T > t + CFL*heff*SetUp.hs){
	    		dt = CFL*heff*SetUp.hs;
	    	}
	    	else {
	    		dt = T - t;
	    	}
			Hermite.interpolate_pd(u,v);
			for(int i = 1; i <= nr; i++){
				for(int j = 1; j <= ns; j++){
					Hermite.recursion_curvilinear(Hermite.ud_interp,Hermite.vd_interp,i,j,
						SetUp.a01d,SetUp.a10d,SetUp.a11d,SetUp.a12d,SetUp.a22d,ud,vd,SetUp.hr,SetUp.hs,dt);
				}
			}
			Hermite.interpolate_dp(ud,vd);
			Hermite.boundaryConditions(m,SetUp.hr,SetUp.hs,ud,vd,SetUp.a01,SetUp.a10,
				SetUp.a11,SetUp.a12,SetUp.a22);
			Hermite.BoundaryPeriodic(ud,vd);
			for(int i = 0; i <= nr; i++){
				for(int j = 0; j <= ns; j++){
					Hermite.recursion_curvilinear(Hermite.u_interp,Hermite.v_interp,i,j,
						SetUp.a01,SetUp.a10,SetUp.a11,SetUp.a12,SetUp.a22,u,v,SetUp.hr,SetUp.hs,dt);
				}
			}
			t += dt;
		 //    Hermite.interpolate_pd(u,v); // Interpolate primal to dual
			// Hermite.oversample(Hermite.ud_interp,U_refined,nsample_r,nsample_s); // Oversample solution	
	  //   	sprintf(fName, "uplot%05d.ext",plotCOUNT);
   //  		Hermite.outPutMat(U_refined,0,nsample_r*nr,0,nsample_s*ns,fName);
			cout << "time = " << t << " Max = " << u.getMax() << endl;
			// plotCOUNT += 1;
	    }
	    cout << endl;
	    Hermite.interpolate_pd(u,v); // Interpolate primal to dual
		Hermite.oversample(Hermite.ud_interp,U_refined,nsample_r,nsample_s); // Oversample solution
		Hermite.oversample(SetUp.Xd,X_refined,nsample_r,nsample_s);
		Hermite.oversample(SetUp.Yd,Y_refined,nsample_r,nsample_s);
    	sprintf(fName, "urefined%05d.ext",COUNT);
    	Hermite.outPutMat(U_refined,0,nsample_r*nr,0,nsample_s*ns,fName);
    	sprintf(fName, "Xrefined%05d.ext",COUNT);
    	Hermite.outPutMat(X_refined,0,nsample_r*nr,0,nsample_s*ns,fName);
    	sprintf(fName, "Yrefined%05d.ext",COUNT);
    	Hermite.outPutMat(Y_refined,0,nsample_r*nr,0,nsample_s*ns,fName);
		COUNT += 1;
	}

	return 0;
}
