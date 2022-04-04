#include "SetUp.h"

void SetUp::compute_grids(){
    r.define(0,nr);
    s.define(0,ns);
    hr = 1.0/double(nr);
    hs = 1.0/double(ns);
    for(int i = 0; i <= nr; i++){
        r(i) = 0.0 + i*hr;
    }
    for(int j = 0; j <= ns; j++){
        s(j) = 0.0 + j*hs;
    }
    // Dual Grids
    rd.define(1,nr);
    sd.define(1,ns);
    for(int i = 1; i <= nr; i++){
        rd(i) = 0.0 - 0.5*hr + hr*double(i);
    }
    for(int j = 1; j <= ns; j++){
        sd(j) = 0.0 - 0.5*hs + hs*double(j);
    }
    // Dual Grid extended
    rd_ext.define(0,nr+1);
    sd_ext.define(0,ns+1);
    for(int i = 0; i <= nr+1; i++){
        rd_ext(i) = 0.0 - 0.5*hr + hr*double(i);
    }
    for(int j = 0; j <= ns+1; j++){
        sd_ext(j) = 0.0 - 0.5*hs + hs*double(j);
    }
}

double cBessel(double nu,double x){
    return boost::math::cyl_bessel_j(nu,x);
}

void SetUp::compute_initial_data(Darray4 &u, Darray4 &v){
    Darray2 U_loc;
    U_loc.define(0,m,0,m);
    U_loc.set_value(0.0);
    double *U_loc_ptr = U_loc.c_ptr();
    int I = 0;
    double r_loc,s_loc;
    for(int j = 0; j<= ns; j++){
        for(int i = 0;i<=nr; i++){
            r_loc = r(i);
            s_loc = s(j);
            indat_(U_loc_ptr,&r_loc,&s_loc,&hr,&hs,&m,&I);
            for(int ids = 0; ids <= m; ids ++){
                for(int idr = 0; idr <= m; idr ++){
                    u(idr,ids,i,j) = U_loc(idr,ids);
                }
            }
        }
    }
}



void transformation_x(int m,int nr_start,int ns_start,int nr, int ns, 
    double hr, double hs,Darray1 &r,Darray1 &s,Darray4 &X){
    Darray2 X_loc;
    int m2p5 = 2*m+5;
    X_loc.define(0,m2p5,0,m2p5);
    X_loc.set_value(0.0);
    double *X_loc_ptr = X_loc.c_ptr();
    int I = 1;
    double r_loc,s_loc;
    for(int j = ns_start; j<= ns; j++){
        for(int i = nr_start;i<=nr; i++){
            r_loc = r(i);
            s_loc = s(j);
            indat_(X_loc_ptr,&r_loc,&s_loc,&hr,&hs,&m2p5,&I);
            for(int ids = 0; ids <= m2p5; ids ++){
                for(int idr = 0; idr <= m2p5; idr ++){
                    X(idr,ids,i,j) = X_loc(idr,ids);
                }
            }
        }
    }
}

void transformation_y(int m,int nr_start,int ns_start,int nr, int ns, 
    double hr, double hs,Darray1 &r,Darray1 &s,Darray4 &Y){
    Darray2 Y_loc;
    int m2p5 = 2*m+5;
    Y_loc.define(0,m2p5,0,m2p5);
    Y_loc.set_value(0.0);
    double *Y_loc_ptr = Y_loc.c_ptr();
    int I = 2;
    double r_loc,s_loc;
    for(int j = ns_start; j<= ns; j++){
        for(int i = nr_start;i<=nr; i++){
            r_loc = r(i);
            s_loc = s(j);
            indat_(Y_loc_ptr,&r_loc,&s_loc,&hr,&hs,&m2p5,&I);
            for(int ids = 0; ids <= m2p5; ids ++){
                for(int idr = 0; idr <= m2p5; idr ++){
                    Y(idr,ids,i,j) = Y_loc(idr,ids);
                }
            }
        }
    }
}

void SetUp::transformation(){
    int m2p5 = 2*m+5;
    // Primal Grid
    X.define(0,m2p5,0,m2p5,0,nr,0,ns);
    Y.define(0,m2p5,0,m2p5,0,nr,0,ns);
    Y.set_value(0.0);
    X.set_value(0.0);
    transformation_x(m,0,0,nr,ns,hr,hs,r,s,X);
    transformation_y(m,0,0,nr,ns,hr,hs,r,s,Y);
    // Dual Grid
    Xd.define(0,m2p5,0,m2p5,1,nr,1,ns);
    Yd.define(0,m2p5,0,m2p5,1,nr,1,ns);
    transformation_x(m,1,1,nr,ns,hr,hs,rd,sd,Xd);
    transformation_y(m,1,1,nr,ns,hr,hs,rd,sd,Yd);
}

void metricInfo(Darray1 &r, Darray1 &s,int m,
    double hr,double hs,int nr_start,int ns_start,
    int nr, int ns,Darray4& X,Darray4 &Y,
    Darray4 &rx,Darray4 &ry,Darray4 &sx, Darray4 &sy){
    Darray2 xr,xs,yr,ys;
    int m2p3 = 2*m+3;
    xr.define(0,m2p3,0,m2p3);
    xs.define(0,m2p3,0,m2p3);
    yr.define(0,m2p3,0,m2p3);
    ys.define(0,m2p3,0,m2p3);
    xr.set_value(0.0);
    xs.set_value(0.0);
    yr.set_value(0.0);
    ys.set_value(0.0);
    double *xr_ptr = xr.c_ptr();
    double *xs_ptr = xs.c_ptr();
    double *yr_ptr = yr.c_ptr();
    double *ys_ptr = ys.c_ptr();
    Darray2 xrys, yrxs;
    xrys.define(0,m2p3,0,m2p3);
    yrxs.define(0,m2p3,0,m2p3);
    xrys.set_value(0.0);
    yrxs.set_value(0.0);
    double *xrys_ptr = xrys.c_ptr();
    double *yrxs_ptr = yrxs.c_ptr();
    for(int i = nr_start; i <= nr; i++){
        for(int j = ns_start; j <= ns; j++){
            // Derivative wrt r
            for(int idr = 0; idr <= m2p3; idr++){
                for(int ids = 0; ids <= m2p3; ids++){
                    xr(idr,ids) = double(idr+1)*X(idr+1,ids,i,j)/hr;
                    yr(idr,ids) = double(idr+1)*Y(idr+1,ids,i,j)/hr;
                }
            }
            // Derivative wrt s
            for(int idr = 0; idr <= m2p3; idr++){
                for(int ids = 0; ids <= m2p3; ids++){
                    xs(idr,ids) = double(ids+1)*X(idr,ids+1,i,j)/hs;
                    ys(idr,ids) = double(ids+1)*Y(idr,ids+1,i,j)/hs;
                }
            }
            
            tpolymul2_(xr_ptr,ys_ptr,xrys_ptr,&m2p3);
            tpolymul2_(yr_ptr,xs_ptr,yrxs_ptr,&m2p3);
            Darray2 J;
            J.define(0,m2p3,0,m2p3);
            J.set_value(0.0);
            double *J_ptr = J.c_ptr();
            for(int idr = 0; idr <= m2p3; idr++){
                for(int ids = 0; ids <= m2p3; ids ++){
                    J(idr,ids) = xrys(idr,ids) - yrxs(idr,ids);
                }
            }   

            Darray2 Denom;
            Denom.define(0,m2p3,0,m2p3);
            Denom.set_value(0.0);
            double polypower;
            polypower = -1.0;
            double *Denom_ptr = Denom.c_ptr();
            tpolypow2_(J_ptr,Denom_ptr,&polypower,&m2p3);
            
            Darray2 Quotient;
            Quotient.define(0,m2p3,0,m2p3);
            Quotient.set_value(0.0);
            double* Quotient_ptr = Quotient.c_ptr();
   
            tpolymul2_(Denom_ptr,ys_ptr,Quotient_ptr,&m2p3);
            for(int idr = 0; idr <= m2p3; idr++){
                for(int ids = 0; ids <= m2p3; ids++){
                    rx(idr,ids,i,j) = Quotient(idr,ids);
                }
            }
            tpolymul2_(Denom_ptr,xs_ptr,Quotient_ptr,&m2p3);
            for(int idr = 0; idr <= m2p3; idr++){
                for(int ids = 0; ids <= m2p3; ids++){
                    ry(idr,ids,i,j) = -Quotient(idr,ids);
                }
            }

            tpolymul2_(Denom_ptr,yr_ptr,Quotient_ptr,&m2p3);
            for(int idr = 0; idr <= m2p3; idr++){
                for(int ids = 0; ids <= m2p3; ids++){
                    sx(idr,ids,i,j) = -Quotient(idr,ids);
                }
            }

            tpolymul2_(Denom_ptr,xr_ptr,Quotient_ptr,&m2p3);
            for(int idr = 0; idr <= m2p3; idr++){
                for(int ids = 0; ids <= m2p3; ids++){
                    sy(idr,ids,i,j) = Quotient(idr,ids);
                }
            }
        }
    }
}


void SetUp::Metric_Info(){
    int m2p3 = 2*m+3;
    rx.define(0,m2p3,0,m2p3,0,nr,0,ns);
    ry.define(0,m2p3,0,m2p3,0,nr,0,ns);
    sx.define(0,m2p3,0,m2p3,0,nr,0,ns);
    sy.define(0,m2p3,0,m2p3,0,nr,0,ns);
    rx.set_value(0.0);
    ry.set_value(0.0);
    sx.set_value(0.0);
    sy.set_value(0.0);
    metricInfo(r,s,m,hr,hs,0,0,nr,ns,X,Y,rx,ry,sx,sy);
    rxd.define(0,m2p3,0,m2p3,1,nr,1,ns);
    ryd.define(0,m2p3,0,m2p3,1,nr,1,ns);
    sxd.define(0,m2p3,0,m2p3,1,nr,1,ns);
    syd.define(0,m2p3,0,m2p3,1,nr,1,ns);
    rxd.set_value(0.0);
    ryd.set_value(0.0);
    sxd.set_value(0.0);
    syd.set_value(0.0);
    metricInfo(rd,sd,m,hr,hs,1,1,nr,ns,Xd,Yd,rxd,ryd,sxd,syd);
}

void computeMetricCoefficients(int i, int j,Darray2 &a01,Darray2 &a10,Darray2 &a11,Darray2 &a12,Darray2 &a22,
    double hr, double hs,Darray4 &Rx, Darray4 &Ry,Darray4 & Sx,Darray4 &Sy,int &m){
    // Set metric info
    int m2p3 = 2*m+3;
    Darray2 rx,ry,sx,sy;
    rx.define(0,m2p3,0,m2p3);
    ry.define(0,m2p3,0,m2p3);
    sx.define(0,m2p3,0,m2p3);
    sy.define(0,m2p3,0,m2p3);
    rx.set_value(0.0);
    ry.set_value(0.0);
    sx.set_value(0.0);
    sy.set_value(0.0);
    // Define pointers to metric data
    double *rx_ptr = rx.c_ptr();
    double *ry_ptr = ry.c_ptr();
    double *sx_ptr = sx.c_ptr();
    double *sy_ptr = sy.c_ptr();

    // Define arrays needed for computation
    Darray2 rx_rx,ry_ry,rx_sx,ry_sy,sx_sx,sy_sy;
    rx_rx.define(0,m2p3,0,m2p3);
    ry_ry.define(0,m2p3,0,m2p3);
    rx_sx.define(0,m2p3,0,m2p3);
    ry_sy.define(0,m2p3,0,m2p3);
    sx_sx.define(0,m2p3,0,m2p3);
    sy_sy.define(0,m2p3,0,m2p3);

    rx_rx.set_value(0.0);
    ry_ry.set_value(0.0);
    rx_sx.set_value(0.0);
    ry_sy.set_value(0.0);
    sx_sx.set_value(0.0);
    sy_sy.set_value(0.0);

    // Define pointers to arrays above
    double *rx_rx_ptr = rx_rx.c_ptr();
    double *ry_ry_ptr = ry_ry.c_ptr();
    double *rx_sx_ptr = rx_sx.c_ptr();
    double *ry_sy_ptr = ry_sy.c_ptr();
    double *sx_sx_ptr = sx_sx.c_ptr();
    double *sy_sy_ptr = sy_sy.c_ptr();

    for(int k=0;k<=m2p3;k++){
        for(int l=0;l<=m2p3;l++){
            rx(k,l) = Rx(k,l,i,j);
        }
    }
    for(int k=0;k<=m2p3;k++){
        for(int l=0;l<=m2p3;l++){
            ry(k,l) = Ry(k,l,i,j);
        }
    }
    
    for(int k=0;k<=m2p3;k++){
        for(int l=0;l<=m2p3;l++){
            sx(k,l) = Sx(k,l,i,j);
        }
    }
    for(int k=0;k<=m2p3;k++){
        for(int l=0;l<=m2p3;l++){
            sy(k,l) = Sy(k,l,i,j);
        }
    }

    tpolymul2_(rx_ptr,rx_ptr,rx_rx_ptr,&m2p3);
    tpolymul2_(ry_ptr,ry_ptr,ry_ry_ptr,&m2p3);
    tpolymul2_(rx_ptr,sx_ptr,rx_sx_ptr,&m2p3);
    tpolymul2_(ry_ptr,sy_ptr,ry_sy_ptr,&m2p3);
    tpolymul2_(sx_ptr,sx_ptr,sx_sx_ptr,&m2p3);
    tpolymul2_(sy_ptr,sy_ptr,sy_sy_ptr,&m2p3);


    // Second Derivatives Compute rxx
    Darray2 rxr, rxr_rx, rxs, rxs_sx;
    rxr.define(0,m2p3,0,m2p3);
    rxr_rx.define(0,m2p3,0,m2p3);
    rxs.define(0,m2p3,0,m2p3);
    rxs_sx.define(0,m2p3,0,m2p3);
    rxr.set_value(0.0);
    rxr_rx.set_value(0.0);
    rxs.set_value(0.0);
    rxs_sx.set_value(0.0);

    // Pointers
    double* rxr_ptr = rxr.c_ptr();
    double* rxr_rx_ptr = rxr_rx.c_ptr();
    double* rxs_ptr = rxs.c_ptr();
    double* rxs_sx_ptr = rxs_sx.c_ptr();

    // Compute rxr, rxs
    for(int idr = 0; idr <= m2p3-1; idr++){
        for(int ids = 0; ids <= m2p3-2; ids++){
            rxr(idr,ids) = (idr+1)*rx(idr+1,ids)/hr;
        }
    }
    for(int idr = 0; idr <= m2p3-1; idr++){
        for(int ids = 0; ids <= m2p3-2; ids++){
            rxs(idr,ids) = (ids+1)*rx(idr,ids+1)/hs;
        }
    }

    tpolymul2_(rxr_ptr,rx_ptr,rxr_rx_ptr,&m2p3);
    tpolymul2_(rxs_ptr,sx_ptr,rxs_sx_ptr,&m2p3);
    Darray2 rxx;
    rxx.define(0,m2p3,0,m2p3);
    rxx.set_value(0.0);
    for(int idr = 0; idr <= m2p3; idr++){
        for(int ids = 0; ids <= m2p3; ids++){
            rxx(idr,ids) = rxr_rx(idr,ids) + rxs_sx(idr,ids);
        }
    }
    // Second Derivatives Compute ryy
    Darray2 ryr, ryr_ry, rys, rys_sy;
    ryr.define(0,m2p3,0,m2p3);
    ryr_ry.define(0,m2p3,0,m2p3);
    rys.define(0,m2p3,0,m2p3);
    rys_sy.define(0,m2p3,0,m2p3);
    ryr.set_value(0.0);
    ryr_ry.set_value(0.0);
    rys.set_value(0.0);
    rys_sy.set_value(0.0);

    // Pointers
    double* ryr_ptr = ryr.c_ptr();
    double* ryr_ry_ptr = ryr_ry.c_ptr();
    double* rys_ptr = rys.c_ptr();
    double* rys_sy_ptr = rys_sy.c_ptr();

    // Compute ryr, rys
    for(int idr = 0; idr <= m2p3-2; idr++){
        for(int ids = 0; ids <= m2p3-1; ids++){
            ryr(idr,ids) = double(idr+1)*ry(idr+1,ids)/hr;
        }
    }
    for(int idr = 0; idr <= m2p3-1; idr++){
        for(int ids = 0; ids <= m2p3-2; ids++){
            rys(idr,ids) = double(ids+1)*ry(idr,ids+1)/hs;
        }
    }

    tpolymul2_(ryr_ptr,ry_ptr,ryr_ry_ptr,&m2p3);
    tpolymul2_(rys_ptr,sy_ptr,rys_sy_ptr,&m2p3);

    
    Darray2 ryy;
    ryy.define(0,m2p3,0,m2p3);
    ryy.set_value(0.0);
    for(int idr = 0; idr <= m2p3; idr++){
        for(int ids = 0; ids <= m2p3; ids++){
            ryy(idr,ids) = ryr_ry(idr,ids) + rys_sy(idr,ids);
        }
    }

    // Second Derivatives Compute sxx
    Darray2 sxr, sxr_rx, sxs, sxs_sx;
    sxr.define(0,m2p3,0,m2p3);
    sxr_rx.define(0,m2p3,0,m2p3);
    sxs.define(0,m2p3,0,m2p3);
    sxs_sx.define(0,m2p3,0,m2p3);
    sxr.set_value(0.0);
    sxr_rx.set_value(0.0);
    sxs.set_value(0.0);
    sxs_sx.set_value(0.0);

    // Pointers
    double* sxr_ptr = sxr.c_ptr();
    double* sxr_rx_ptr = sxr_rx.c_ptr();
    double* sxs_ptr = sxs.c_ptr();
    double* sxs_sx_ptr = sxs_sx.c_ptr();

    // Compute rxr, rxs
    for(int idr = 0; idr <= m2p3-2; idr++){
        for(int ids = 0; ids <= m2p3-1; ids++){
            sxr(idr,ids) = double(idr+1)*sx(idr+1,ids)/hr;
        }
    }
    for(int idr = 0; idr <= m2p3-1; idr++){
        for(int ids = 0; ids <= m2p3-2; ids++){
            sxs(idr,ids) = double(ids+1)*sx(idr,ids+1)/hs;
        }
    }

    tpolymul2_(sxr_ptr,rx_ptr,sxr_rx_ptr,&m2p3);
    tpolymul2_(sxs_ptr,sx_ptr,sxs_sx_ptr,&m2p3);
    
    Darray2 sxx;
    sxx.define(0,m2p3,0,m2p3);
    sxx.set_value(0.0);
    for(int idr = 0; idr <= m2p3; idr++){
        for(int ids = 0; ids <= m2p3; ids++){
            sxx(idr,ids) = sxr_rx(idr,ids) + sxs_sx(idr,ids);
        }
    }

    // Second Derivatives Compute syy
    Darray2 syr, syr_ry, sys, sys_sy;
    syr.define(0,m2p3,0,m2p3);
    syr_ry.define(0,m2p3,0,m2p3);
    sys.define(0,m2p3,0,m2p3);
    sys_sy.define(0,m2p3,0,m2p3);
    syr.set_value(0.0);
    syr_ry.set_value(0.0);
    sys.set_value(0.0);
    sys_sy.set_value(0.0);

    // Pointers
    double* syr_ptr = syr.c_ptr();
    double* syr_ry_ptr = syr_ry.c_ptr();
    double* sys_ptr = sys.c_ptr();
    double* sys_sy_ptr = sys_sy.c_ptr();

    // Compute syr, sys
    for(int idr = 0; idr <= m2p3-2; idr++){
        for(int ids = 0; ids <= m2p3-1; ids++){
            syr(idr,ids) = double(idr+1)*sy(idr+1,ids)/hr;
        }
    }
    for(int idr = 0; idr <= m2p3-1; idr++){
        for(int ids = 0; ids <= m2p3-2; ids++){
            sys(idr,ids) = double(ids+1)*sy(idr,ids+1)/hs;
        }
    }

    tpolymul2_(syr_ptr,ry_ptr,syr_ry_ptr,&m2p3);
    tpolymul2_(sys_ptr,sy_ptr,sys_sy_ptr,&m2p3);
    
    Darray2 syy;
    syy.define(0,m2p3,0,m2p3);
    syy.set_value(0.0);
    for(int idr = 0; idr <= m2p3; idr++){
        for(int ids = 0; ids <= m2p3; ids++){
            syy(idr,ids) = syr_ry(idr,ids) + sys_sy(idr,ids);
        }
    }
    a01.set_value(0.0);
    a10.set_value(0.0);
    a11.set_value(0.0);
    a12.set_value(0.0);
    a22.set_value(0.0);
    
    // Set variable coefficient values
    for(int idr = 0; idr <= m2p3; idr++){
        for(int ids = 0; ids <= m2p3; ids++){
            a01(idr,ids) = sxx(idr,ids) + syy(idr,ids);
            a10(idr,ids) = rxx(idr,ids) + ryy(idr,ids);
            a11(idr,ids) = rx_rx(idr,ids) + ry_ry(idr,ids);
            a12(idr,ids) = 2*(rx_sx(idr,ids) + ry_sy(idr,ids));
            a22(idr,ids) = sx_sx(idr,ids) + sy_sy(idr,ids);
        }
    }
}

void SetUp::metricCoeffs(){
    // Metric Coefficients
    int m2p3 = 2*m+3;
    a01.define(0,m2p3,0,m2p3,0,nr,0,ns);
    a10.define(0,m2p3,0,m2p3,0,nr,0,ns);
    a11.define(0,m2p3,0,m2p3,0,nr,0,ns);
    a12.define(0,m2p3,0,m2p3,0,nr,0,ns);
    a22.define(0,m2p3,0,m2p3,0,nr,0,ns);
    // Metric Coefficients on dual
    a01d.define(0,m2p3,0,m2p3,1,nr,1,ns);
    a10d.define(0,m2p3,0,m2p3,1,nr,1,ns);
    a11d.define(0,m2p3,0,m2p3,1,nr,1,ns);
    a12d.define(0,m2p3,0,m2p3,1,nr,1,ns);
    a22d.define(0,m2p3,0,m2p3,1,nr,1,ns);
    // Metric local
    Darray2 a01_loc, a10_loc, a11_loc, a12_loc, a22_loc;
    a01_loc.define(0,m2p3,0,m2p3);
    a10_loc.define(0,m2p3,0,m2p3);
    a11_loc.define(0,m2p3,0,m2p3);
    a12_loc.define(0,m2p3,0,m2p3);
    a22_loc.define(0,m2p3,0,m2p3);

    for(int j = 0; j <= ns; j++){
        for(int i = 0; i <= nr; i++){
            computeMetricCoefficients(i,j,a01_loc,a10_loc,a11_loc,a12_loc,a22_loc,
                hr,hs,rx,ry,sx,sy,m);
            for(int ids = 0; ids <= m2p3; ids++){
                for(int idr = 0; idr <= m2p3; idr++){
                    a01(idr,ids,i,j) = a01_loc(idr,ids);
                    a10(idr,ids,i,j) = a10_loc(idr,ids);
                    a11(idr,ids,i,j) = a11_loc(idr,ids);
                    a12(idr,ids,i,j) = a12_loc(idr,ids);
                    a22(idr,ids,i,j) = a22_loc(idr,ids);
                }
            }
        }
    }

    for(int j = 1; j <= ns; j++){
        for(int i = 1; i <= nr; i++){
            computeMetricCoefficients(i,j,a01_loc,a10_loc,a11_loc,a12_loc,a22_loc,
                hr,hs,rxd,ryd,sxd,syd,m);
            for(int ids = 0; ids <= m2p3; ids++){
                for(int idr = 0; idr <= m2p3; idr++){
                    a01d(idr,ids,i,j) = a01_loc(idr,ids);
                    a10d(idr,ids,i,j) = a10_loc(idr,ids);
                    a11d(idr,ids,i,j) = a11_loc(idr,ids);
                    a12d(idr,ids,i,j) = a12_loc(idr,ids);
                    a22d(idr,ids,i,j) = a22_loc(idr,ids);
                }
            }
        }
    }
}
// -----------------------------------------------------------------------
double SetUp::computeMinJacobian(double hr,double hs,int nr, int ns,Darray4& X,Darray4 &Y){
    double Jmin = 10e10;
    double J = 10e10;
    double xr,yr,xs,ys;
    int m2p5 = 2*m+5;
    for(int i = 0; i <= nr; i++){
        for(int j = 0; j <= ns; j++){
            // Derivative wrt r
            for(int idr = 0; idr <= m2p5 - 1; idr++){
                for(int ids = 0; ids <= m2p5; ids++){
                    xr = X(1,0,i,j)/hr;
                    yr = Y(1,0,i,j)/hr;
                }
            }
            // Derivative wrt s
            for(int idr = 0; idr <= m2p5; idr++){
                for(int ids = 0; ids <= m2p5 - 1; ids++){
                    xs = X(0,1,i,j)/hs;
                    ys = Y(0,1,i,j)/hs;
                }
            }
            J = pow(abs(xr*ys - yr*xs),0.5);
            Jmin = (J < Jmin) ? J : Jmin;
        }
    }
    return Jmin;
}
// -----------------------------------------------------------------------
