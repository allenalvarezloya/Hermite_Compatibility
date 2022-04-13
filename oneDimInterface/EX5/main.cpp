#include <iostream>
#include <iomanip>
#include <cmath>
#include "Darrays.h"
#include "SetUp.h"
#include "Hermite.h"

using namespace std;

int main(){
	char fName[100];
	for(int k = 1; k <= 5; k++){
		for(int m = 1; m <= 3; m++){
			double CFL = 0.5;
			int n1 = 10+2*k;
			double C = 1.0; // Wave Speed 
			// Allocate memory for solution arrays (Primal)
			Darray2 u1,um1;
			u1.define(0,m,0,n1);
			um1.define(0,m,0,n1);
			u1.set_value(0.0);
			um1.set_value(0.0);
			// Allocate memory for solution arrays (Dual)
			Darray2 ud1,umd1;
			ud1.define(0,m,1,n1);
			umd1.define(0,m-1,1,n1);
			ud1.set_value(0.0);
			umd1.set_value(0.0);
			// Set up the problem
			SetUp SetUp1;
			SetUp1.nr = n1;
			SetUp1.m = m;
			SetUp1.compute_grids();
			double h = SetUp1.hr;
			double T = 1.3;
			double t = 0.0;
			int nsteps;
		    double dt = CFL*h;
		    nsteps = ceil(T/dt);
		    dt = T/nsteps;
			SetUp1.compute_initial_data(u1,SetUp1.r,0,n1,0.0);
			SetUp1.compute_initial_data(ud1,SetUp1.rd,1,n1,-0.5*dt);
			Hermite Hermite1;
			Hermite1.nr = n1;
			Hermite1.m = m;
			Hermite1.compute_Hmaps();
		    for(int steps = 1; steps <= nsteps; steps++){
				Hermite1.interpPD(u1);
				for(int i = 1; i <= n1; i++){				
					Hermite1.recursion2(Hermite1.ud_interp,ud1,i,h,C,dt);			
				}
				Hermite1.interpDP(ud1);
				Hermite1.boundaryConditions(m,SetUp1.hr,ud1,t+0.5*dt);
				for(int i = 0; i <= n1; i++){
					Hermite1.recursion2(Hermite1.u_interp,u1,i,h,C,dt);
				}
				Darray1 S1;
				S1.define(0,n1);
				for(int i = 0; i <= n1; i++){
					S1(i) = u1(0,i);
				}
			    sprintf(fName, "solution1%05d.ext",steps);
			    Hermite1.outPutVec(S1,0,n1,fName);
			    t += dt;
			}
			sprintf(fName,"u1m%doversampled%05d.ext",m,k);
			Hermite1.interpPD(u1);
			Hermite1.oversample(Hermite1.ud_interp,SetUp1.nr,m,fName);
		}
	}
	return 0;
}
