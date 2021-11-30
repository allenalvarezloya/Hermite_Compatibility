#include <iostream>
#include <cmath>
#include "SetUp.h"
#include "Hermite.h"
using namespace std;
extern "C" void dgetrs_(char *,int *, int *, double *, int *, int *, double *, int *, int *);
int main(){
	int m = 2;
	int COUNT = 1;
	int refine = 3;
	Darray1 errors1, errors2, errorsinf;
	errors1.define(1,refine);
	errors2.define(1,refine);
	errorsinf.define(1,refine);
	errors1.set_value(0.0);
	errors2.set_value(0.0);
	errorsinf.set_value(0.0); 
	for(int n = 10; n <= 10*pow(2,refine-1); n *= 2){
		int nr = n;
		int ns = n;
		// // Allocate memory for solution arrays (Primal)
		Darray4 u,v;
		u.define(0,m,0,m,0,nr,0,ns);
		v.define(0,m-1,0,m-1,0,nr,0,ns);
		u.set_value(0.0);
		v.set_value(0.0);
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
		SetUp.transformation_noScale();
		SetUp.compute_initial_data(u,v);
		SetUp.Metric_Info();
		SetUp.metricCoeffs();
		SetUp.Metric_Info_noScale();
		SetUp.metricCoeffs_noScale();

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
		double t, dt, T;
	    T = 2.0/pow(2.0,0.5);
	    t = 0.0;
	    dt = 0.0;
	    while (t < T){
	    	if (T > t + 0.5*SetUp.hs){
	    		dt = 0.5*SetUp.hs;
	    	}
	    	else {
	    		dt = T - t;
	    	}
			Hermite.interpolate_pd(u,v);
			for(int i = 1; i <= nr; i++){
				for(int j = 1; j <= ns; j++){
					Hermite.cartesian_recursion(Hermite.ud_interp,Hermite.vd_interp,ud,vd,i,j,SetUp.hr,SetUp.hs,dt);
				}
			}
			Hermite.interpolate_dp(ud,vd);
			Hermite.boundaryConditions(m,SetUp.hr,SetUp.hs,ud,vd,SetUp.a01_ns,SetUp.a10_ns,
				SetUp.a11_ns,SetUp.a12_ns,SetUp.a22_ns);
			Hermite.BoundaryPeriodic(ud,vd);
			for(int i = 0; i <= nr; i++){
				for(int j = 0; j <= ns; j++){
					Hermite.cartesian_recursion(Hermite.u_interp,Hermite.v_interp,u,v,i,j,SetUp.hr,SetUp.hs,dt);
				}
			}
			t += dt;
			cout << "time = " << t << endl;
	    }
	    Hermite.interpolate_pd(u,v); // Interpolate primal to dual
		Hermite.oversample(Hermite.ud_interp,t,errors1,errors2,errorsinf,COUNT); // Oversample solution
		COUNT += 1;
	}
	FILE *extFile1 = fopen("L1errors.ext", "w");
    for (int i=1; i<=refine; i++){
        fprintf(extFile1, " %18.10e", errors1(i));
        fprintf(extFile1,"\n");
    }   
    fclose(extFile1);
	FILE *extFile2 = fopen("L2errors.ext", "w");
    for (int i=1; i<=refine; i++){
        fprintf(extFile2, " %18.10e", errors2(i));
        fprintf(extFile2,"\n");
    }   
    fclose(extFile2);
    	FILE *extFileinf = fopen("Linferrors.ext", "w");
    for (int i=1; i<=refine; i++){
        fprintf(extFileinf, " %18.10e", errorsinf(i));
        fprintf(extFileinf,"\n");
    }   
    fclose(extFileinf);

	return 0;
}
