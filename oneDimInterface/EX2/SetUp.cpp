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

void SetUp::compute_initial_data(Darray2 &u,Darray2 &v){
    Darray1 U_loc;
    U_loc.define(0,m);
    double *U_loc_ptr = U_loc.c_ptr();
    double r_loc;
    int I;
    I = 0;
    for(int i = 0; i <= nr; i++){
        r_loc = r(i);
        indat_(U_loc_ptr,&r_loc,&hr,&m,&I);
        for(int idr = 0; idr <= m; idr++){
            u(idr,i) = U_loc(idr);
        }
    }
    int m_minus = m-1;
    Darray1 V_loc;
    V_loc.define(0,m-1);
    double *V_loc_ptr = V_loc.c_ptr();
    I = 1;
    for(int i = 0; i <= nr; i++){
        r_loc = r(i);
        indat_(V_loc_ptr,&r_loc,&hr,&m_minus,&I);
        for(int idr = 0; idr <= m_minus; idr++){
            v(idr,i) = V_loc(idr);
        }
    }
}