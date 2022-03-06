#include "SetUp.h"
void SetUp::compute_grids(){
    // Primal Grid
    r.define(0,nr);
    hr = 1.0/double(nr);
    for(int i = 0; i <= nr; i++){
        r(i) = 0.0 + i*hr;
    }
    // Dual Grids
    rd.define(0,nr+1);
    for(int i = 0; i <= nr+1; i++){
        rd(i) = 0.0 - 0.5*hr + hr*double(i);
    }
}

void SetUp::compute_initial_data(Darray2 &u,int I){
    Darray1 U_loc;
    U_loc.define(0,m);
    double *U_loc_ptr = U_loc.c_ptr();
    double r_loc;
    for(int i = 0; i <= nr; i++){
        r_loc = r(i);
        indat_(U_loc_ptr,&r_loc,&hr,&m,&I);
        for(int idr = 0; idr <= m; idr++){
            u(idr,i) = U_loc(idr);
        }
    }
}