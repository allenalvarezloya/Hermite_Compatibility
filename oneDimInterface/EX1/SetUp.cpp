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

void SetUp::compute_initial_data(Darray2 &u,Darray1 &x,int ns,int ne,double t){
    Darray1 U_loc;
    U_loc.define(0,m);
    double *U_loc_ptr = U_loc.c_ptr();
    double r_loc;
    for(int i = ns; i <= ne; i++){
        r_loc = x(i);
        indat_(U_loc_ptr,&r_loc,&t,&hr,&m);
        for(int idr = 0; idr <= m; idr++){
            u(idr,i) = U_loc(idr);
        }
    }
    // double scale = 0.0;
    // int mod_m;
    // for(int i = ns; i <= ne; i++){
    //     for(int idr = 0; idr <= m; idr ++){
    //         scale = pow(2*M_PI*hr,idr)/tgamma(idr+1);
    //         mod_m = idr % 4;
    //         if (mod_m == 0){
    //             u(idr,i) = scale*sin(2*M_PI*x(i))*cos(2*M_PI*t);
    //         }
    //         else if(mod_m == 1){
    //             u(idr,i) = scale*cos(2*M_PI*x(i))*cos(2*M_PI*t);
    //         }
    //         else if(mod_m == 2){
    //             u(idr,i) = -scale*sin(2*M_PI*x(i))*cos(2*M_PI*t);
    //         }
    //         else if(mod_m == 3 ){
    //             u(idr,i) = -scale*cos(2*M_PI*x(i))*cos(2*M_PI*t);
    //         }

    //     }
    // }
}