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
	int refine = 0;
	Darray1 errors1, errors2, errorsinf;
	errors1.define(1,refine);
	errors2.define(1,refine);
	errorsinf.define(1,refine);
	errors1.set_value(0.0);
	errors2.set_value(0.0);
	errorsinf.set_value(0.0); 
	char fName[100];
	int base = 2;
	// For refinement
    int nsample_r, nsample_s;
    nsample_r = 2;    // Oversample in r
    nsample_s = 2;    // Oversample in s
	for(int n = base; n <= 50; n ++){
	    double k31, k33,gamma,alpha;
	    k31 = 6.3801618959239383506237;
	    k33 = 13.01520072169843441983;
	    gamma = k31/k33;
	    alpha = 1-gamma;
		int nr = n;
		int ns = n;


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

		Hermite.computeSingular(m,SetUp.hr,SetUp.hs,SetUp.a01,SetUp.a10,SetUp.a11,SetUp.a12,SetUp.a22,n);

	}

	return 0;
}
