#include <iostream>
#include <cmath>
#include "Darrays.h"
#include "SetUp.h"
#include "Hermite.h"

using namespace std;

int main(){
	char fName[100];
	// for(int k = 1; k <= 5; k++){
		int k = 3;
		int m = 3;
		int n1 = 10*pow(2.0,k-1);
		int n2 = 20*pow(2.0,k-1);
		int nh1 = n1/2;
		int nh2 = n2/2;
		double C = 1.0; // Wave Speed on left
		double D = 2.0; // Wave Speed on right
		// Allocate memory for solution arrays (Primal)
		Darray2 u1,v1;
		u1.define(0,m,0,n1);
		v1.define(0,m-1,0,n1);
		u1.set_value(0.0);
		v1.set_value(0.0);
		// Allocate memory for solution arrays (Dual)
		Darray2 ud1,vd1;
		ud1.define(0,m,1,n1);
		vd1.define(0,m-1,1,n1);
		ud1.set_value(0.0);
		vd1.set_value(0.0);
		// Set up the problem
		int I = 1;
		SetUp SetUp1;
		SetUp1.nr = n1;
		SetUp1.m = m;
		SetUp1.compute_grids();
		SetUp1.compute_initial_data(u1,I);

		Hermite Hermite1;
		Hermite1.nr = n1;
		Hermite1.m = m;
		Hermite1.compute_Hmaps();


		// Second grid! 
			// Allocate memory for solution arrays (Primal)
		Darray2 u2,v2;
		u2.define(0,m,0,n2);
		v2.define(0,m-1,0,n2);
		u2.set_value(0.0);
		v2.set_value(0.0);
		// Allocate memory for solution arrays (Dual)
		Darray2 ud2,vd2;
		ud2.define(0,m,1,n2);
		vd2.define(0,m-1,1,n2);
		ud2.set_value(0.0);
		vd2.set_value(0.0);
		// Set up the problem
		I = 0;
		SetUp SetUp2;
		SetUp2.nr = n2;
		SetUp2.m = m;
		SetUp2.compute_grids();
		SetUp2.compute_initial_data(u2,I);


		Hermite Hermite2;
		Hermite2.nr = n2;
		Hermite2.m = m;
		Hermite2.compute_Hmaps();
		for(int i = nh1+1; i <= n1; i++){
			for(int idr = 0; idr <= m; idr++){
				u1(idr,i) = 0.0;
			}
		}
		for(int i = nh1+1; i <= n1; i++){
			for(int idr = 0; idr <= m-1; idr++){
				v1(idr,i) = 0.0;
			}
		}

		for(int i = 0; i <= nh2 - 1; i++){
			for(int idr = 0; idr <= m; idr++){
				u2(idr,i) = 0.0;
			}
		}
		for(int i = 0; i <= nh2 - 1; i++){
			for(int idr = 0; idr <= m-1; idr++){
				v2(idr,i) = 0.0;
			}
		}
		// end Second grid

		int COUNT = 0;
		double time, dt, T;
	    T = 2.0;
	    time = 0.0;
	    while (time < T){
	    	if (T > time + SetUp2.hr/D){
	    		dt = SetUp2.hr/D;
	    	}
	    	else {
	    		dt = T - time;
	    	}
			Hermite1.interpPD(u1,v1);
			Hermite2.interpPD(u2,v2);
			for(int i = 1; i <= nh1; i++){
				Hermite1.recursion(Hermite1.ud_interp,Hermite1.vd_interp,ud1,vd1,i,SetUp1.hr,C,dt);
			}
			for(int i = nh2+1; i <= n2; i++){
				Hermite2.recursion(Hermite2.ud_interp,Hermite2.vd_interp,ud2,vd2,i,SetUp2.hr,D,dt);
			}
			Hermite1.interpDP(ud1,vd1);
			Hermite2.interpDP(ud2,vd2);
			Hermite1.boundaryConditions(m,SetUp1.hr,ud1,vd1);
			Hermite2.boundaryConditions(m,SetUp2.hr,ud2,vd2);
			Hermite1.interfaceConditions(m,SetUp1.hr,SetUp2.hr,C,D,ud1,ud2,vd1,vd2,nh1,nh2+1,nh1,1);
			Hermite2.interfaceConditions(m,SetUp1.hr,SetUp2.hr,C,D,ud1,ud2,vd1,vd2,nh1,nh2+1,nh2,2);
			for(int i = 0; i <= nh1; i++){
				Hermite1.recursion(Hermite1.u_interp,Hermite1.v_interp,u1,v1,i,SetUp1.hr,C,dt);
			}
			for(int i = nh2; i <= n2; i++){
				Hermite2.recursion(Hermite2.u_interp,Hermite2.v_interp,u2,v2,i,SetUp2.hr,D,dt);
			}
			time += dt;
			Darray1 S1;
			S1.define(0,n1);
			for(int i = 0; i <= n1; i++){
				S1(i) = u1(0,i);
			}
		    sprintf(fName, "solution1%05d.ext",COUNT);
		    Hermite1.outPutVec(S1,0,n1,fName);
			Darray1 S2;
			S2.define(0,n2);
			for(int i = 0; i <= n2; i++){
				S2(i) = u2(0,i);
			}
		    sprintf(fName, "solution2%05d.ext",COUNT);
		    Hermite1.outPutVec(S2,0,n2,fName);
		    COUNT += 1;
		}
		// sprintf(fName,"u1oversampled%05d.ext",k);
		// Hermite1.interpPD(u1,v1);
		// Hermite1.oversample(Hermite1.ud_interp,SetUp1.nr,m,fName);
		// sprintf(fName,"u2oversampled%05d.ext",k);
		// Hermite2.interpPD(u2,v2);
		// Hermite2.oversample(Hermite2.ud_interp,SetUp2.nr,m,fName);
	// }

	return 0;
}
