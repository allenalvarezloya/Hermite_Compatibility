#include <iostream>
#include <cmath>
#include "Darrays.h"
#include "SetUp.h"
#include "Hermite.h"

using namespace std;

int main(){
	char fName[100];
	for(int k = 2; k <= 10; k++){
		for(int m = 1; m <=3; m++){
			int n = 10+4*k;
			double C = 1.0; // Wave Speed 
			// Allocate memory for solution arrays (Primal)
			Darray2 u,v;
			u.define(0,m,0,n);
			v.define(0,m-1,0,n);
			u.set_value(0.0);
			v.set_value(0.0);
			// Allocate memory for solution arrays (Dual)
			Darray2 ud,vd;
			ud.define(0,m,1,n);
			vd.define(0,m-1,1,n);
			ud.set_value(0.0);
			vd.set_value(0.0);
			// Set up the problem
			SetUp SetUp;
			SetUp.nr = n;
			SetUp.m = m;
			SetUp.compute_grids();
			SetUp.compute_initial_data(u,v);

			Hermite Hermite;
			Hermite.nr = n;
			Hermite.m = m;
			Hermite.compute_Hmaps();


			int COUNT = 0;
			double CFL = 0.5;
			double h = SetUp.hr;
		    double T = 1.3;
			int nsteps;
		    double dt = CFL*h;
		    nsteps = ceil(T/dt);
		    dt = T/nsteps;
		    for(int steps = 1; steps <= nsteps; steps++){
				Hermite.interpPD(u,v);
				for(int i = 1; i <= n; i++){
					Hermite.recursion(Hermite.ud_interp,Hermite.vd_interp,ud,vd,i,h,C,dt);
				}
				Hermite.interpDP(ud,vd);
				Hermite.boundaryConditions(m,SetUp.hr,ud,vd);
				for(int i = 0; i <= n; i++){
					Hermite.recursion(Hermite.u_interp,Hermite.v_interp,u,v,i,h,C,dt);
				}
				// Darray1 S;
				// S.define(0,n);
				// for(int i = 0; i <= n; i++){
				// 	S(i) = u(0,i);
				// }
			 //    sprintf(fName, "solution%05d.ext",COUNT);
			 //    Hermite.outPutVec(S,0,n,fName);
			 //    COUNT += 1;
			}
			sprintf(fName,"um%doversampled%05d.ext",m,k);
			Hermite.interpPD(u,v);
			Hermite.oversample(Hermite.ud_interp,SetUp.nr,m,fName);
		}
	}
	return 0;
}
