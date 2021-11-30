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
    double k31, k33,gamma,alpha,rpart;
    rpart = 1.0;
    k31 = 6.3801618959239383506237;
    k33 = 13.01520072169843441983;
    gamma = k31/k33;
    alpha = 1-gamma;
    double scale;
    int case_r,case_s;
    for(int i = 0; i <= nr; i++){
        for(int j = 0; j <= ns; j++){
            for(int idr = 0; idr <= m; idr++){
                scale = pow(hr,idr)/tgamma(idr+1);
                case_r = idr;
                if(case_r == 0){
                    rpart = cBessel(3,(alpha*r(i)+gamma)*k33);
                }
                else if(case_r == 1){
                    rpart = scale*0.5*alpha*k33*(cBessel(2,(alpha*r(i)+gamma)*k33)-cBessel(4,(alpha*r(i)+gamma)*k33));
                }
                else if(case_r == 2){
                    rpart = scale*0.5*alpha*k33*(0.5*alpha*k33*(cBessel(1,(alpha*r(i)+gamma)*k33)-cBessel(3,(alpha*r(i)+gamma)*k33))
                        -0.5*alpha*k33*(cBessel(3,(alpha*r(i)+gamma)*k33)-cBessel(5,(alpha*r(i)+gamma)*k33)));
                }
                else if(case_r == 3){
                    rpart = scale*pow(0.5*k33*alpha,3)*(cBessel(0,(alpha*r(i)+gamma)*k33)-3*cBessel(2,(alpha*r(i)+gamma)*k33)
                        +3*cBessel(4,(alpha*r(i)+gamma)*k33)-cBessel(6,(alpha*r(i)+gamma)*k33));
                }
                for(int ids = 0; ids <= m; ids++){
                    scale = pow(3*2*M_PI*hs,ids)/tgamma(ids+1);
                    case_s = ids % 4;
                    if(case_s == 0){
                        u(idr,ids,i,j) = scale*rpart*cos(3*2*M_PI*s(j));
                    }
                    else if(case_s == 1){
                        u(idr,ids,i,j) = -scale*rpart*sin(3*2*M_PI*s(j));
                    }
                    else if(case_s == 2){
                        u(idr,ids,i,j) = -scale*rpart*cos(3*2*M_PI*s(j));
                    }
                    else if(case_s == 3){
                        u(idr,ids,i,j) = scale*rpart*sin(3*2*M_PI*s(j));
                    }
                }
            }
        }
    }
}

void SetUp::compute_initial_data_dual_ext(Darray4 &u, Darray4 &v){
    double k31, k33,gamma,alpha,rpart;
    rpart = 1.0;
    k31 = 6.3801618959239383506237;
    k33 = 13.01520072169843441983;
    gamma = k31/k33;
    alpha = 1-gamma;
    double scale;
    int case_r,case_s;
    for(int i = 0; i <= nr+1; i++){
        for(int j = 0; j <= ns+1; j++){
            for(int idr = 0; idr <= m; idr++){
                scale = pow(hr,idr)/tgamma(idr+1);
                case_r = idr;
                if(case_r == 0){
                    rpart = cBessel(3,(alpha*rd_ext(i)+gamma)*k33);
                }
                else if(case_r == 1){
                    rpart = scale*0.5*alpha*k33*(cBessel(2,(alpha*rd_ext(i)+gamma)*k33)-cBessel(4,(alpha*rd_ext(i)+gamma)*k33));
                }
                else if(case_r == 2){
                    rpart = scale*0.5*alpha*k33*(0.5*alpha*k33*(cBessel(1,(alpha*rd_ext(i)+gamma)*k33)-cBessel(3,(alpha*rd_ext(i)+gamma)*k33))
                        -0.5*alpha*k33*(cBessel(3,(alpha*rd_ext(i)+gamma)*k33)-cBessel(5,(alpha*rd_ext(i)+gamma)*k33)));
                }
                else if(case_r == 3){
                    rpart = scale*pow(0.5*k33*alpha,3)*(cBessel(0,(alpha*rd_ext(i)+gamma)*k33)-3*cBessel(2,(alpha*rd_ext(i)+gamma)*k33)
                        +3*cBessel(4,(alpha*rd_ext(i)+gamma)*k33)-cBessel(6,(alpha*rd_ext(i)+gamma)*k33));
                }
                for(int ids = 0; ids <= m; ids++){
                    scale = pow(3*2*M_PI*hs,ids)/tgamma(ids+1);
                    case_s = ids % 4;
                    if(case_s == 0){
                        u(idr,ids,i,j) = scale*rpart*cos(3*2*M_PI*sd_ext(j));
                    }
                    else if(case_s == 1){
                        u(idr,ids,i,j) = -scale*rpart*sin(3*2*M_PI*sd_ext(j));
                    }
                    else if(case_s == 2){
                        u(idr,ids,i,j) = -scale*rpart*cos(3*2*M_PI*sd_ext(j));
                    }
                    else if(case_s == 3){
                        u(idr,ids,i,j) = scale*rpart*sin(3*2*M_PI*sd_ext(j));
                    }
                }
            }
        }
    }
}

void transformation_x(int m,int nr_start,int ns_start,int nr, int ns, 
	double hr, double hs,Darray1 &r,Darray1 &s,Darray4 &X,int metricOrder){
    int case_s;
    double scale;
    double k31, k33,gamma,alpha;
    k31 = 6.3801618959239383506237;
    k33 = 13.01520072169843441983;
    gamma = k31/k33;
    alpha = 1-gamma;
    for(int j = ns_start; j<= ns; j++){
        for(int i = nr_start;i<=nr; i++){
            for(int idr=0;idr<=metricOrder;idr++){
                for(int ids=0;ids<=metricOrder;ids++){
                    scale = (pow(alpha*hr,idr)*pow(2*M_PI*hs,ids))/(tgamma(idr+1)*tgamma(ids+1));
                    if(idr == 0){
                        case_s = ids % 4;
                        if(case_s == 0){
                            X(idr,ids,i,j) = scale*(alpha*r(i)+gamma)*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 1){
                            X(idr,ids,i,j) = -scale*(alpha*r(i)+gamma)*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 2){
                            X(idr,ids,i,j) = -scale*(alpha*r(i)+gamma)*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 3){
                            X(idr,ids,i,j) = scale*(alpha*r(i)+gamma)*sin(2*M_PI*s(j));
                        }
                    }
                    else if(idr == 1){
                        case_s = ids % 4;
                        if(case_s == 0){
                            X(idr,ids,i,j) = scale*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 1){
                            X(idr,ids,i,j) = -scale*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 2){
                            X(idr,ids,i,j) = -scale*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 3){
                            X(idr,ids,i,j) = scale*sin(2*M_PI*s(j));
                        }
                    }
                }
            }
        }
    }
}

void transformation_y(int m,int nr_start,int ns_start,int nr, int ns, 
	double hr, double hs,Darray1 &r,Darray1 &s,Darray4 &Y,int metricOrder){
    int case_s;
    double scale;
    double k31, k33,gamma,alpha;
    k31 = 6.3801618959239383506237;
    k33 = 13.01520072169843441983;
    gamma = k31/k33;
    alpha = 1-gamma;
    for(int j = ns_start; j<= ns; j++){
        for(int i = nr_start;i<=nr; i++){
            for(int idr=0;idr<=metricOrder;idr++){
                for(int ids=0;ids<=metricOrder;ids++){
                    scale = (pow(alpha*hr,idr)*pow(2*M_PI*hs,ids))/(tgamma(idr+1)*tgamma(ids+1));
                    if(idr == 0){
                        case_s = ids % 4;
                        if(case_s == 0){
                            Y(idr,ids,i,j) = scale*(alpha*r(i)+gamma)*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 1){
                            Y(idr,ids,i,j) = scale*(alpha*r(i)+gamma)*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 2){
                            Y(idr,ids,i,j) = -scale*(alpha*r(i)+gamma)*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 3){
                            Y(idr,ids,i,j) = -scale*(alpha*r(i)+gamma)*cos(2*M_PI*s(j));
                        }
                    }
                    else if(idr == 1){
                        case_s = ids % 4;
                        if(case_s == 0){
                            Y(idr,ids,i,j) = scale*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 1){
                            Y(idr,ids,i,j) = scale*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 2){
                            Y(idr,ids,i,j) = -scale*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 3){
                            Y(idr,ids,i,j) = -scale*cos(2*M_PI*s(j));
                        }
                    }
                }
            }
        }
    }
}

void SetUp::transformation(){
    // Primal Grid
    X.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    Y.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    Y.set_value(0.0);
    X.set_value(0.0);
    transformation_x(m,0,0,nr,ns,hr,hs,r,s,X,metricOrder);
    transformation_y(m,0,0,nr,ns,hr,hs,r,s,Y,metricOrder);
    // Dual Grid
    Xd.define(0,metricOrder,0,metricOrder,1,nr,1,ns);
    Yd.define(0,metricOrder,0,metricOrder,1,nr,1,ns);
    transformation_x(m,1,1,nr,ns,hr,hs,rd,sd,Xd,metricOrder);
    transformation_y(m,1,1,nr,ns,hr,hs,rd,sd,Yd,metricOrder);
}

void metricInfo(Darray1 &r, Darray1 &s,int m,
    double hr,double hs,int nr_start,int ns_start,
    int nr, int ns,Darray4& X,Darray4 &Y,
    Darray4 &rx,Darray4 &ry,Darray4 &sx, Darray4 &sy,int metricOrder){
    Darray2 xr,xs,yr,ys;
    xr.define(0,metricOrder,0,metricOrder);
    xs.define(0,metricOrder,0,metricOrder);
    yr.define(0,metricOrder,0,metricOrder);
    ys.define(0,metricOrder,0,metricOrder);
    xr.set_value(0.0);
    xs.set_value(0.0);
    yr.set_value(0.0);
    ys.set_value(0.0);
    double *xr_ptr = xr.c_ptr();
    double *xs_ptr = xs.c_ptr();
    double *yr_ptr = yr.c_ptr();
    double *ys_ptr = ys.c_ptr();
    Darray2 xrys, yrxs;
    xrys.define(0,metricOrder,0,metricOrder);
    yrxs.define(0,metricOrder,0,metricOrder);
    xrys.set_value(0.0);
    yrxs.set_value(0.0);
    double *xrys_ptr = xrys.c_ptr();
    double *yrxs_ptr = yrxs.c_ptr();
    for(int i = nr_start; i <= nr; i++){
        for(int j = ns_start; j <= ns; j++){
            // Derivative wrt r
            for(int idr = 0; idr <= metricOrder - 1; idr++){
                for(int ids = 0; ids <= metricOrder; ids++){
                    xr(idr,ids) = double(idr+1)*X(idr+1,ids,i,j)/hr;
                    yr(idr,ids) = double(idr+1)*Y(idr+1,ids,i,j)/hr;
                }
            }
            // Derivative wrt s
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder - 1; ids++){
                    xs(idr,ids) = double(ids+1)*X(idr,ids+1,i,j)/hs;
                    ys(idr,ids) = double(ids+1)*Y(idr,ids+1,i,j)/hs;
                }
            }
            
            tpolymul2_(xr_ptr,ys_ptr,xrys_ptr,&metricOrder);
            tpolymul2_(yr_ptr,xs_ptr,yrxs_ptr,&metricOrder);
            Darray2 J;
            J.define(0,metricOrder,0,metricOrder);
            J.set_value(0.0);
            double *J_ptr = J.c_ptr();
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder; ids ++){
                    J(idr,ids) = xrys(idr,ids) - yrxs(idr,ids);
                }
            }   
            Darray2 Denom;
            Denom.define(0,metricOrder,0,metricOrder);
            Denom.set_value(0.0);
            double polypower;
            polypower = -1.0;
            double *Denom_ptr = Denom.c_ptr();
            tpolypow2_(J_ptr,Denom_ptr,&polypower,&metricOrder);
            
            Darray2 Quotient;
            Quotient.define(0,metricOrder,0,metricOrder);
            Quotient.set_value(0.0);
            double* Quotient_ptr = Quotient.c_ptr();
   
            tpolymul2_(Denom_ptr,ys_ptr,Quotient_ptr,&metricOrder);
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder; ids++){
                    rx(idr,ids,i,j) = Quotient(idr,ids);
                }
            }
            tpolymul2_(Denom_ptr,xs_ptr,Quotient_ptr,&metricOrder);
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder; ids++){
                    ry(idr,ids,i,j) = -Quotient(idr,ids);
                }
            }

            tpolymul2_(Denom_ptr,yr_ptr,Quotient_ptr,&metricOrder);
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder; ids++){
                    sx(idr,ids,i,j) = -Quotient(idr,ids);
                }
            }

            tpolymul2_(Denom_ptr,xr_ptr,Quotient_ptr,&metricOrder);
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder; ids++){
                    sy(idr,ids,i,j) = Quotient(idr,ids);
                }
            }
        }
    }
}

void SetUp::Metric_Info(){
    rx.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    ry.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    sx.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    sy.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    rx.set_value(0.0);
    ry.set_value(0.0);
    sx.set_value(0.0);
    sy.set_value(0.0);
    metricInfo(r,s,m,hr,hs,0,0,nr,ns,X,Y,rx,ry,sx,sy,metricOrder);
    rxd.define(0,metricOrder,0,metricOrder,1,nr,1,ns);
    ryd.define(0,metricOrder,0,metricOrder,1,nr,1,ns);
    sxd.define(0,metricOrder,0,metricOrder,1,nr,1,ns);
    syd.define(0,metricOrder,0,metricOrder,1,nr,1,ns);
    rxd.set_value(0.0);
    ryd.set_value(0.0);
    sxd.set_value(0.0);
    syd.set_value(0.0);
    metricInfo(rd,sd,m,hr,hs,1,1,nr,ns,Xd,Yd,rxd,ryd,sxd,syd,metricOrder);
}

void computeMetricCoefficients(int i, int j,Darray2 &a01,Darray2 &a10,Darray2 &a11,Darray2 &a12,Darray2 &a22,
    int metricOrder,double hr, double hs,Darray4 &Rx, Darray4 &Ry,Darray4 & Sx,Darray4 &Sy){
    // Set metric info
    Darray2 rx,ry,sx,sy;
    rx.define(0,metricOrder,0,metricOrder);
    ry.define(0,metricOrder,0,metricOrder);
    sx.define(0,metricOrder,0,metricOrder);
    sy.define(0,metricOrder,0,metricOrder);
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
    rx_rx.define(0,metricOrder,0,metricOrder);
    ry_ry.define(0,metricOrder,0,metricOrder);
    rx_sx.define(0,metricOrder,0,metricOrder);
    ry_sy.define(0,metricOrder,0,metricOrder);
    sx_sx.define(0,metricOrder,0,metricOrder);
    sy_sy.define(0,metricOrder,0,metricOrder);

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

    for(int k=0;k<=metricOrder;k++){
        for(int l=0;l<=metricOrder;l++){
            rx(k,l) = Rx(k,l,i,j);
        }
    }
    for(int k=0;k<=metricOrder;k++){
        for(int l=0;l<=metricOrder;l++){
            ry(k,l) = Ry(k,l,i,j);
        }
    }
    
    for(int k=0;k<=metricOrder;k++){
        for(int l=0;l<=metricOrder;l++){
            sx(k,l) = Sx(k,l,i,j);
        }
    }
    for(int k=0;k<=metricOrder;k++){
        for(int l=0;l<=metricOrder;l++){
            sy(k,l) = Sy(k,l,i,j);
        }
    }

    tpolymul2_(rx_ptr,rx_ptr,rx_rx_ptr,&metricOrder);
    tpolymul2_(ry_ptr,ry_ptr,ry_ry_ptr,&metricOrder);
    tpolymul2_(rx_ptr,sx_ptr,rx_sx_ptr,&metricOrder);
    tpolymul2_(ry_ptr,sy_ptr,ry_sy_ptr,&metricOrder);
    tpolymul2_(sx_ptr,sx_ptr,sx_sx_ptr,&metricOrder);
    tpolymul2_(sy_ptr,sy_ptr,sy_sy_ptr,&metricOrder);


    // Second Derivatives Compute rxx
    Darray2 rxr, rxr_rx, rxs, rxs_sx;
    rxr.define(0,metricOrder,0,metricOrder);
    rxr_rx.define(0,metricOrder,0,metricOrder);
    rxs.define(0,metricOrder,0,metricOrder);
    rxs_sx.define(0,metricOrder,0,metricOrder);
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
    for(int idr = 0; idr <= metricOrder - 1; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            rxr(idr,ids) = (idr+1)*rx(idr+1,ids)/hr;
        }
    }
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder-1; ids++){
            rxs(idr,ids) = (ids+1)*rx(idr,ids+1)/hs;
        }
    }

    tpolymul2_(rxr_ptr,rx_ptr,rxr_rx_ptr,&metricOrder);
    tpolymul2_(rxs_ptr,sx_ptr,rxs_sx_ptr,&metricOrder);
    Darray2 rxx;
    rxx.define(0,metricOrder,0,metricOrder);
    rxx.set_value(0.0);
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            rxx(idr,ids) = rxr_rx(idr,ids) + rxs_sx(idr,ids);
        }
    }
    // Second Derivatives Compute ryy
    Darray2 ryr, ryr_ry, rys, rys_sy;
    ryr.define(0,metricOrder,0,metricOrder);
    ryr_ry.define(0,metricOrder,0,metricOrder);
    rys.define(0,metricOrder,0,metricOrder);
    rys_sy.define(0,metricOrder,0,metricOrder);
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
    for(int idr = 0; idr <= metricOrder - 1; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            ryr(idr,ids) = double(idr+1)*ry(idr+1,ids)/hr;
        }
    }
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder-1; ids++){
            rys(idr,ids) = double(ids+1)*ry(idr,ids+1)/hs;
        }
    }

    tpolymul2_(ryr_ptr,ry_ptr,ryr_ry_ptr,&metricOrder);
    tpolymul2_(rys_ptr,sy_ptr,rys_sy_ptr,&metricOrder);

    
    Darray2 ryy;
    ryy.define(0,metricOrder,0,metricOrder);
    ryy.set_value(0.0);
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            ryy(idr,ids) = ryr_ry(idr,ids) + rys_sy(idr,ids);
        }
    }

    // Second Derivatives Compute sxx
    Darray2 sxr, sxr_rx, sxs, sxs_sx;
    sxr.define(0,metricOrder,0,metricOrder);
    sxr_rx.define(0,metricOrder,0,metricOrder);
    sxs.define(0,metricOrder,0,metricOrder);
    sxs_sx.define(0,metricOrder,0,metricOrder);
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
    for(int idr = 0; idr <= metricOrder - 1; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            sxr(idr,ids) = double(idr+1)*sx(idr+1,ids)/hr;
        }
    }
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder-1; ids++){
            sxs(idr,ids) = double(ids+1)*sx(idr,ids+1)/hs;
        }
    }

    tpolymul2_(sxr_ptr,rx_ptr,sxr_rx_ptr,&metricOrder);
    tpolymul2_(sxs_ptr,sx_ptr,sxs_sx_ptr,&metricOrder);
    
    Darray2 sxx;
    sxx.define(0,metricOrder,0,metricOrder);
    sxx.set_value(0.0);
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            sxx(idr,ids) = sxr_rx(idr,ids) + sxs_sx(idr,ids);
        }
    }

    // Second Derivatives Compute syy
    Darray2 syr, syr_ry, sys, sys_sy;
    syr.define(0,metricOrder,0,metricOrder);
    syr_ry.define(0,metricOrder,0,metricOrder);
    sys.define(0,metricOrder,0,metricOrder);
    sys_sy.define(0,metricOrder,0,metricOrder);
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
    for(int idr = 0; idr <= metricOrder - 1; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            syr(idr,ids) = double(idr+1)*sy(idr+1,ids)/hr;
        }
    }
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder-1; ids++){
            sys(idr,ids) = double(ids+1)*sy(idr,ids+1)/hs;
        }
    }

    tpolymul2_(syr_ptr,ry_ptr,syr_ry_ptr,&metricOrder);
    tpolymul2_(sys_ptr,sy_ptr,sys_sy_ptr,&metricOrder);
    
    Darray2 syy;
    syy.define(0,metricOrder,0,metricOrder);
    syy.set_value(0.0);
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            syy(idr,ids) = syr_ry(idr,ids) + sys_sy(idr,ids);
        }
    }
    a01.set_value(0.0);
    a10.set_value(0.0);
    a11.set_value(0.0);
    a12.set_value(0.0);
    a22.set_value(0.0);
    
    // Set variable coefficient values
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
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
    a01.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    a10.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    a11.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    a12.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    a22.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    // Metric Coefficients on dual
    a01d.define(0,metricOrder,0,metricOrder,1,nr,1,ns);
    a10d.define(0,metricOrder,0,metricOrder,1,nr,1,ns);
    a11d.define(0,metricOrder,0,metricOrder,1,nr,1,ns);
    a12d.define(0,metricOrder,0,metricOrder,1,nr,1,ns);
    a22d.define(0,metricOrder,0,metricOrder,1,nr,1,ns);
    // Metric local
    Darray2 a01_loc, a10_loc, a11_loc, a12_loc, a22_loc;
    a01_loc.define(0,metricOrder,0,metricOrder);
    a10_loc.define(0,metricOrder,0,metricOrder);
    a11_loc.define(0,metricOrder,0,metricOrder);
    a12_loc.define(0,metricOrder,0,metricOrder);
    a22_loc.define(0,metricOrder,0,metricOrder);

    for(int j = 0; j <= ns; j++){
        for(int i = 0; i <= nr; i++){
            computeMetricCoefficients(i,j,a01_loc,a10_loc,a11_loc,a12_loc,a22_loc,
                metricOrder,hr,hs,rx,ry,sx,sy);
            for(int ids = 0; ids <= metricOrder; ids++){
                for(int idr = 0; idr <= metricOrder; idr++){
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
                metricOrder,hr,hs,rxd,ryd,sxd,syd);
            for(int ids = 0; ids <= metricOrder; ids++){
                for(int idr = 0; idr <= metricOrder; idr++){
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



void transformation_x_noScale(int m,int nr_start,int ns_start,int nr, int ns,
    Darray1 &r,Darray1 &s,Darray4 &X,int metricOrder){
    int case_s;
    double scale;
    double k31, k33,gamma,alpha;
    k31 = 6.3801618959239383506237;
    k33 = 13.01520072169843441983;
    gamma = k31/k33;
    alpha = 1-gamma;
    for(int j = ns_start; j<= ns; j++){
        for(int i = nr_start;i<=nr; i++){
            for(int idr=0;idr<=metricOrder;idr++){
                for(int ids=0;ids<=metricOrder;ids++){
                    scale = (pow(alpha,idr)*pow(2*M_PI,ids))/(tgamma(idr+1)*tgamma(ids+1));
                    if(idr == 0){
                        case_s = ids % 4;
                        if(case_s == 0){
                            X(idr,ids,i,j) = scale*(alpha*r(i)+gamma)*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 1){
                            X(idr,ids,i,j) = -scale*(alpha*r(i)+gamma)*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 2){
                            X(idr,ids,i,j) = -scale*(alpha*r(i)+gamma)*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 3){
                            X(idr,ids,i,j) = scale*(alpha*r(i)+gamma)*sin(2*M_PI*s(j));
                        }
                    }
                    else if(idr == 1){
                        case_s = ids % 4;
                        if(case_s == 0){
                            X(idr,ids,i,j) = scale*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 1){
                            X(idr,ids,i,j) = -scale*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 2){
                            X(idr,ids,i,j) = -scale*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 3){
                            X(idr,ids,i,j) = scale*sin(2*M_PI*s(j));
                        }
                    }
                }
            }
        }
    }
}

void transformation_y_noScale(int m,int nr_start,int ns_start,int nr, int ns, 
    Darray1 &r,Darray1 &s,Darray4 &Y,int metricOrder){
    int case_s;
    double scale;
    double k31, k33,gamma,alpha;
    k31 = 6.3801618959239383506237;
    k33 = 13.01520072169843441983;
    gamma = k31/k33;
    alpha = 1-gamma;
    for(int j = ns_start; j<= ns; j++){
        for(int i = nr_start;i<=nr; i++){
            for(int idr=0;idr<=metricOrder;idr++){
                for(int ids=0;ids<=metricOrder;ids++){
                    scale = (pow(alpha,idr)*pow(2*M_PI,ids))/(tgamma(idr+1)*tgamma(ids+1));
                    if(idr == 0){
                        case_s = ids % 4;
                        if(case_s == 0){
                            Y(idr,ids,i,j) = scale*(alpha*r(i)+gamma)*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 1){
                            Y(idr,ids,i,j) = scale*(alpha*r(i)+gamma)*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 2){
                            Y(idr,ids,i,j) = -scale*(alpha*r(i)+gamma)*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 3){
                            Y(idr,ids,i,j) = -scale*(alpha*r(i)+gamma)*cos(2*M_PI*s(j));
                        }
                    }
                    else if(idr == 1){
                        case_s = ids % 4;
                        if(case_s == 0){
                            Y(idr,ids,i,j) = scale*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 1){
                            Y(idr,ids,i,j) = scale*cos(2*M_PI*s(j));
                        }
                        else if(case_s == 2){
                            Y(idr,ids,i,j) = -scale*sin(2*M_PI*s(j));
                        }
                        else if(case_s == 3){
                            Y(idr,ids,i,j) = -scale*cos(2*M_PI*s(j));
                        }
                    }
                }
            }
        }
    }
}

void SetUp::transformation_noScale(){
    // Primal Grid
    X_ns.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    Y_ns.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    Y_ns.set_value(0.0);
    X_ns.set_value(0.0);
    transformation_x_noScale(m,0,0,nr,ns,r,s,X_ns,metricOrder);
    transformation_y_noScale(m,0,0,nr,ns,r,s,Y_ns,metricOrder);
}

void metricInfo_noScale(Darray1 &r, Darray1 &s,int m,
    int nr_start,int ns_start,
    int nr, int ns,Darray4& X,Darray4 &Y,
    Darray4 &rx,Darray4 &ry,Darray4 &sx, Darray4 &sy,int metricOrder){
    Darray2 xr,xs,yr,ys;
    xr.define(0,metricOrder,0,metricOrder);
    xs.define(0,metricOrder,0,metricOrder);
    yr.define(0,metricOrder,0,metricOrder);
    ys.define(0,metricOrder,0,metricOrder);
    xr.set_value(0.0);
    xs.set_value(0.0);
    yr.set_value(0.0);
    ys.set_value(0.0);
    double *xr_ptr = xr.c_ptr();
    double *xs_ptr = xs.c_ptr();
    double *yr_ptr = yr.c_ptr();
    double *ys_ptr = ys.c_ptr();
    Darray2 xrys, yrxs;
    xrys.define(0,metricOrder,0,metricOrder);
    yrxs.define(0,metricOrder,0,metricOrder);
    xrys.set_value(0.0);
    yrxs.set_value(0.0);
    double *xrys_ptr = xrys.c_ptr();
    double *yrxs_ptr = yrxs.c_ptr();
    for(int i = nr_start; i <= nr; i++){
        for(int j = ns_start; j <= ns; j++){
            // Derivative wrt r
            for(int idr = 0; idr <= metricOrder - 1; idr++){
                for(int ids = 0; ids <= metricOrder; ids++){
                    xr(idr,ids) = double(idr+1)*X(idr+1,ids,i,j);
                    yr(idr,ids) = double(idr+1)*Y(idr+1,ids,i,j);
                }
            }
            // Derivative wrt s
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder - 1; ids++){
                    xs(idr,ids) = double(ids+1)*X(idr,ids+1,i,j);
                    ys(idr,ids) = double(ids+1)*Y(idr,ids+1,i,j);
                }
            }
            
            tpolymul2_(xr_ptr,ys_ptr,xrys_ptr,&metricOrder);
            tpolymul2_(yr_ptr,xs_ptr,yrxs_ptr,&metricOrder);
            Darray2 J;
            J.define(0,metricOrder,0,metricOrder);
            J.set_value(0.0);
            double *J_ptr = J.c_ptr();
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder; ids ++){
                    J(idr,ids) = xrys(idr,ids) - yrxs(idr,ids);
                }
            }   
            Darray2 Denom;
            Denom.define(0,metricOrder,0,metricOrder);
            Denom.set_value(0.0);
            double polypower;
            polypower = -1.0;
            double *Denom_ptr = Denom.c_ptr();
            tpolypow2_(J_ptr,Denom_ptr,&polypower,&metricOrder);
            
            Darray2 Quotient;
            Quotient.define(0,metricOrder,0,metricOrder);
            Quotient.set_value(0.0);
            double* Quotient_ptr = Quotient.c_ptr();
   
            tpolymul2_(Denom_ptr,ys_ptr,Quotient_ptr,&metricOrder);
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder; ids++){
                    rx(idr,ids,i,j) = Quotient(idr,ids);
                }
            }
            tpolymul2_(Denom_ptr,xs_ptr,Quotient_ptr,&metricOrder);
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder; ids++){
                    ry(idr,ids,i,j) = -Quotient(idr,ids);
                }
            }

            tpolymul2_(Denom_ptr,yr_ptr,Quotient_ptr,&metricOrder);
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder; ids++){
                    sx(idr,ids,i,j) = -Quotient(idr,ids);
                }
            }

            tpolymul2_(Denom_ptr,xr_ptr,Quotient_ptr,&metricOrder);
            for(int idr = 0; idr <= metricOrder; idr++){
                for(int ids = 0; ids <= metricOrder; ids++){
                    sy(idr,ids,i,j) = Quotient(idr,ids);
                }
            }
        }
    }
}

void SetUp::Metric_Info_noScale(){
    rx_ns.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    ry_ns.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    sx_ns.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    sy_ns.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    rx_ns.set_value(0.0);
    ry_ns.set_value(0.0);
    sx_ns.set_value(0.0);
    sy_ns.set_value(0.0);
    metricInfo_noScale(r,s,m,0,0,nr,ns,X_ns,Y_ns,rx_ns,ry_ns,sx_ns,sy_ns,metricOrder);

}

void computeMetricCoefficients_noScale(int i, int j,Darray2 &a01,Darray2 &a10,Darray2 &a11,Darray2 &a12,Darray2 &a22,
    int metricOrder,Darray4 &Rx, Darray4 &Ry,Darray4 & Sx,Darray4 &Sy){
    // Set metric info
    Darray2 rx,ry,sx,sy;
    rx.define(0,metricOrder,0,metricOrder);
    ry.define(0,metricOrder,0,metricOrder);
    sx.define(0,metricOrder,0,metricOrder);
    sy.define(0,metricOrder,0,metricOrder);
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
    rx_rx.define(0,metricOrder,0,metricOrder);
    ry_ry.define(0,metricOrder,0,metricOrder);
    rx_sx.define(0,metricOrder,0,metricOrder);
    ry_sy.define(0,metricOrder,0,metricOrder);
    sx_sx.define(0,metricOrder,0,metricOrder);
    sy_sy.define(0,metricOrder,0,metricOrder);

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

    for(int k=0;k<=metricOrder;k++){
        for(int l=0;l<=metricOrder;l++){
            rx(k,l) = Rx(k,l,i,j);
        }
    }
    for(int k=0;k<=metricOrder;k++){
        for(int l=0;l<=metricOrder;l++){
            ry(k,l) = Ry(k,l,i,j);
        }
    }
    
    for(int k=0;k<=metricOrder;k++){
        for(int l=0;l<=metricOrder;l++){
            sx(k,l) = Sx(k,l,i,j);
        }
    }
    for(int k=0;k<=metricOrder;k++){
        for(int l=0;l<=metricOrder;l++){
            sy(k,l) = Sy(k,l,i,j);
        }
    }

    tpolymul2_(rx_ptr,rx_ptr,rx_rx_ptr,&metricOrder);
    tpolymul2_(ry_ptr,ry_ptr,ry_ry_ptr,&metricOrder);
    tpolymul2_(rx_ptr,sx_ptr,rx_sx_ptr,&metricOrder);
    tpolymul2_(ry_ptr,sy_ptr,ry_sy_ptr,&metricOrder);
    tpolymul2_(sx_ptr,sx_ptr,sx_sx_ptr,&metricOrder);
    tpolymul2_(sy_ptr,sy_ptr,sy_sy_ptr,&metricOrder);


    // Second Derivatives Compute rxx
    Darray2 rxr, rxr_rx, rxs, rxs_sx;
    rxr.define(0,metricOrder,0,metricOrder);
    rxr_rx.define(0,metricOrder,0,metricOrder);
    rxs.define(0,metricOrder,0,metricOrder);
    rxs_sx.define(0,metricOrder,0,metricOrder);
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
    for(int idr = 0; idr <= metricOrder - 1; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            rxr(idr,ids) = (idr+1)*rx(idr+1,ids);
        }
    }
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder-1; ids++){
            rxs(idr,ids) = (ids+1)*rx(idr,ids+1);
        }
    }

    tpolymul2_(rxr_ptr,rx_ptr,rxr_rx_ptr,&metricOrder);
    tpolymul2_(rxs_ptr,sx_ptr,rxs_sx_ptr,&metricOrder);
    Darray2 rxx;
    rxx.define(0,metricOrder,0,metricOrder);
    rxx.set_value(0.0);
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            rxx(idr,ids) = rxr_rx(idr,ids) + rxs_sx(idr,ids);
        }
    }
    // Second Derivatives Compute ryy
    Darray2 ryr, ryr_ry, rys, rys_sy;
    ryr.define(0,metricOrder,0,metricOrder);
    ryr_ry.define(0,metricOrder,0,metricOrder);
    rys.define(0,metricOrder,0,metricOrder);
    rys_sy.define(0,metricOrder,0,metricOrder);
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
    for(int idr = 0; idr <= metricOrder - 1; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            ryr(idr,ids) = double(idr+1)*ry(idr+1,ids);
        }
    }
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder-1; ids++){
            rys(idr,ids) = double(ids+1)*ry(idr,ids+1);
        }
    }

    tpolymul2_(ryr_ptr,ry_ptr,ryr_ry_ptr,&metricOrder);
    tpolymul2_(rys_ptr,sy_ptr,rys_sy_ptr,&metricOrder);

    
    Darray2 ryy;
    ryy.define(0,metricOrder,0,metricOrder);
    ryy.set_value(0.0);
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            ryy(idr,ids) = ryr_ry(idr,ids) + rys_sy(idr,ids);
        }
    }

    // Second Derivatives Compute sxx
    Darray2 sxr, sxr_rx, sxs, sxs_sx;
    sxr.define(0,metricOrder,0,metricOrder);
    sxr_rx.define(0,metricOrder,0,metricOrder);
    sxs.define(0,metricOrder,0,metricOrder);
    sxs_sx.define(0,metricOrder,0,metricOrder);
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
    for(int idr = 0; idr <= metricOrder - 1; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            sxr(idr,ids) = double(idr+1)*sx(idr+1,ids);
        }
    }
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder-1; ids++){
            sxs(idr,ids) = double(ids+1)*sx(idr,ids+1);
        }
    }

    tpolymul2_(sxr_ptr,rx_ptr,sxr_rx_ptr,&metricOrder);
    tpolymul2_(sxs_ptr,sx_ptr,sxs_sx_ptr,&metricOrder);
    
    Darray2 sxx;
    sxx.define(0,metricOrder,0,metricOrder);
    sxx.set_value(0.0);
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            sxx(idr,ids) = sxr_rx(idr,ids) + sxs_sx(idr,ids);
        }
    }

    // Second Derivatives Compute syy
    Darray2 syr, syr_ry, sys, sys_sy;
    syr.define(0,metricOrder,0,metricOrder);
    syr_ry.define(0,metricOrder,0,metricOrder);
    sys.define(0,metricOrder,0,metricOrder);
    sys_sy.define(0,metricOrder,0,metricOrder);
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
    for(int idr = 0; idr <= metricOrder - 1; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            syr(idr,ids) = double(idr+1)*sy(idr+1,ids);
        }
    }
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder-1; ids++){
            sys(idr,ids) = double(ids+1)*sy(idr,ids+1);
        }
    }

    tpolymul2_(syr_ptr,ry_ptr,syr_ry_ptr,&metricOrder);
    tpolymul2_(sys_ptr,sy_ptr,sys_sy_ptr,&metricOrder);
    
    Darray2 syy;
    syy.define(0,metricOrder,0,metricOrder);
    syy.set_value(0.0);
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            syy(idr,ids) = syr_ry(idr,ids) + sys_sy(idr,ids);
        }
    }
    a01.set_value(0.0);
    a10.set_value(0.0);
    a11.set_value(0.0);
    a12.set_value(0.0);
    a22.set_value(0.0);
    
    // Set variable coefficient values
    for(int idr = 0; idr <= metricOrder; idr++){
        for(int ids = 0; ids <= metricOrder; ids++){
            a01(idr,ids) = sxx(idr,ids) + syy(idr,ids);
            a10(idr,ids) = rxx(idr,ids) + ryy(idr,ids);
            a11(idr,ids) = rx_rx(idr,ids) + ry_ry(idr,ids);
            a12(idr,ids) = 2*(rx_sx(idr,ids) + ry_sy(idr,ids));
            a22(idr,ids) = sx_sx(idr,ids) + sy_sy(idr,ids);
        }
    }
}

void SetUp::metricCoeffs_noScale(){
    // Metric Coefficients
    a01_ns.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    a10_ns.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    a11_ns.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    a12_ns.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    a22_ns.define(0,metricOrder,0,metricOrder,0,nr,0,ns);
    // Metric local
    Darray2 a01_loc, a10_loc, a11_loc, a12_loc, a22_loc;
    a01_loc.define(0,metricOrder,0,metricOrder);
    a10_loc.define(0,metricOrder,0,metricOrder);
    a11_loc.define(0,metricOrder,0,metricOrder);
    a12_loc.define(0,metricOrder,0,metricOrder);
    a22_loc.define(0,metricOrder,0,metricOrder);

    for(int j = 0; j <= ns; j++){
        for(int i = 0; i <= nr; i++){
            computeMetricCoefficients_noScale(i,j,a01_loc,a10_loc,a11_loc,a12_loc,a22_loc,
                metricOrder,rx_ns,ry_ns,sx_ns,sy_ns);
            for(int ids = 0; ids <= metricOrder; ids++){
                for(int idr = 0; idr <= metricOrder; idr++){
                    a01_ns(idr,ids,i,j) = a01_loc(idr,ids);
                    a10_ns(idr,ids,i,j) = a10_loc(idr,ids);
                    a11_ns(idr,ids,i,j) = a11_loc(idr,ids);
                    a12_ns(idr,ids,i,j) = a12_loc(idr,ids);
                    a22_ns(idr,ids,i,j) = a22_loc(idr,ids);
                }
            }
        }
    }

}



