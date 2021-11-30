#include <iostream>
#include <cmath>
#include "SetUp.h"
#include "Hermite.h"
using namespace std;
extern "C" void dgetrs_(char *,int *, int *, double *, int *, int *, double *, int *, int *);
int main(){
	int m = 3;
	int n = 80;
	int nr = n;
	int ns = 4*n;
	// Allocate memory for solution arrays (Primal)
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
	// Extended
	Darray4 ud_ext, vd_ext;
	ud_ext.define(0,m,0,m,0,nr+1,0,ns+1);
	vd_ext.define(0,m-1,0,m-1,0,nr+1,0,ns+1);
	ud_ext.set_value(0.0);
	vd_ext.set_value(0.0);
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
	SetUp.compute_initial_data_dual_ext(ud_ext,vd_ext);
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
    double k31, k33,gamma,alpha;
    k31 = 6.3801618959239383506237;
    k33 = 13.01520072169843441983;
    gamma = k31/k33;
    alpha = 1-gamma;
	double t, dt, T;
    T = 0.5*SetUp.hs;
    t = 0.0;
    dt = 0.0;
    int COUNT = 1;
    while (t < T){
    	if (T > t + 0.5*SetUp.hs){
    		dt = 0.5*SetUp.hs;
    	}
    	else {
    		dt = T - t;
    	}
		Hermite.interpolate_pd(u,v);
		Hermite.interpolate_dp_ext(ud_ext,vd_ext);
		Hermite.boundaryConditions(m,SetUp.hr,SetUp.hs,Hermite.ud_interp,Hermite.vd_interp,SetUp.a01_ns,SetUp.a10_ns,
			SetUp.a11_ns,SetUp.a12_ns,SetUp.a22_ns);
		for(int i = 0; i <= nr; i++){
			for(int j = 0; j <= ns; j++){
				Hermite.recursion_curvilinear(Hermite.u_interp,Hermite.v_interp,i,j,
					SetUp.a01,SetUp.a10,SetUp.a11,SetUp.a12,SetUp.a22,u,v,SetUp.hr,SetUp.hs,dt);
			}
		}
		t += dt;
	    // Hermite.interpolate_pd(u,v); // Interpolate primal to dual
    	// Hermite.oversample(Hermite.ud_interp,t,COUNT); // Oversample solution
    	// COUNT += 1;
    }
    Hermite.interpolate_pd(u,v); // Interpolate primal to dual
	Hermite.oversample(Hermite.ud_interp,0.5*t,COUNT); // Oversample solution
	return 0;
}
