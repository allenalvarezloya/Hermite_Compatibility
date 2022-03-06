#include "Hermite.h"
using namespace std;
//--------------------------------------------------
void Hermite::compute_Hmaps() {
    // This function allocates space for Hmat_u and Hmat_v
    // It then calls on the HermiteMap routine to fill in
    // the entries
    Hmat_u.define(0,2*m+1,0,2*m+1);
    Hmat_v.define(0,2*m-1,0,2*m-1);
    Hmat_u.set_value(0.0);
    Hmat_v.set_value(0.0);
    HermiteMap(m,0,1,0.5,0,Hmat_u);
    HermiteMap(m-1,0,1,0.5,0,Hmat_v);
}

void Hermite::binomial(int m,Darray2 &coeffs){
    //#% ! Computes the binomial coefficients of order up through m
    //#% IMPLICIT NONE
    //#% INTEGER, INTENT(IN) :: m
    //#% DOUBLE PRECISION, DIMENSION(0:m,0:m), INTENT(OUT) :: coeffs
    //#% INTEGER :: i,j
    coeffs(0,0) = 1.0;
    
    for(int i=1;i<=m;i++){
        coeffs(i,0) = 1.0;
        for(int j=1;j<=i;j++){
            coeffs(i,j) = coeffs(i,j-1)*(i-j+1)/j;
        }
    }
}
//--------------------------------------------------
void Hermite::HermiteMap(int m,double xl,double xr,double xc,int icase,Darray2 &tmat){



    //#% ! This subroutine computes the coefficient matrix tmat which
    //#% ! transfers derivative data of order through m at xl and xr
    //#% ! to the function and derivative data at xc
    //#% ! That is, if p is a polynomial of degree 2m+1
    //#% !
    //#% !  h^k D^k p (xc)/k! = sum_(j=0)^m tmat(k,j) h^j D^j p(xl)/j!
    //#% !                                + tmat(k,j+m+1) h^j D^j p(xr)/j!
    //#% !
    //#% !  icase < 0 => xc=xl (left boundary case)
    //#% !  icase = 0 => xl < xc < xr
    //#% !  icase > 0 => xc=xr (right boundary case)
    //#% IMPLICIT NONE
    //#% INTEGER, INTENT(IN) :: m,icase
    //#% DOUBLE PRECISION, INTENT(IN) :: xl,xr,xc
    //#% DOUBLE PRECISION, DIMENSION(0:2*m+1,0:2*m+1), INTENT(OUT) :: tmat
    //#% DOUBLE PRECISION :: h,z,zc,adl,adr,sign,c1l,c1r,c2l,c2r
    //#% INTEGER :: i,j,k
    //
    //#% ! Compute in normalized coordinates

    double h, z, zc;
    h = xr-xl;
    z = (xc-xl)/h;
    zc = z-1.0;
    
    if (icase > 0){
        z = 1.0;
        zc = 0.0;
    }
    else if (icase < 0){
        z=0.0;
        zc=-1.0;
    }
    Darray2 bcofs;
    bcofs.define(0,m+1,0,m+1);
    bcofs.set_value(0.0);
    int mp = m+1;
    binomial(mp,bcofs);

    
    //    #% !
    //    #% ! We begin with the Hermite-Lagrange interpolants:
    //    #% !
    //    #% !   Q_j (z) = z^(m+1) sum_{k=j}^m h_{kj} (z-1)^k,
    //    #% !
    //    #% !   j=0, ... , m
    //    #% !
    //    #% !   satisfying Q_j = (z-1)^j + O((z-1)^(m+1)) at z=1
    //    #% !
    //    #% !   After some algebra once can show:
    //    #% !
    //    #% !   h_jj = 1,  h_kj = -sum_{p=j}^{k-1} b_(k-p)^(m+1) h_pj ,
    //    #% !              for k>j
    //    #% !
    //    #% !   here b_(k-p)^(m+1) is the binomial coefficient (m+1)!/((k-p)!(m+1-k+p)!)
    //    #% !
    //    #% ! To construct the matrix we
    //    #% ! now evaluate the interpolant and its derivatives at z
    //    #% !
    //    #% ! Data to the left is handled by a simple change of variables z -> 1-z
    //    #% !
    //    #% ! Start with the last column - note that we directly use the recursive
    //    #% ! definition of the h's to fold previously computed columns into old
    //    #% ! ones. Note that the polynomial is now centered about the midpoint
    //    #% !
    

    for(int i = 0; i<=2*m+1;i++){
        double adl, adr, c1l, c1r, c2l, c2r;
        adl = 0.0;
        adr = 0.0;
        for(int j = std::max(0,i-m);j<=std::min(i,m+1);j++){
            if((m-i+j) == 0){
                c2l = 1.0;
                c2r = 1.0;
            }
            else{
                c2l = pow(z,m-i+j);
                c2r = pow(zc,m-i+j);
            }
            if((m+1-j) == 0){
                c1l = 1.0;
                c1r = 1.0;
            }
            else{
                c1l = pow(zc,m+1-j);
                c1r = pow(z,m+1-j);
            }
            adr += bcofs(m+1,j)*bcofs(m,i-j)*c1r*c2r;
            adl += bcofs(m+1,j)*bcofs(m,i-j)*c1l*c2l;
        }
        tmat(i,2*m+1)=adr;
        tmat(i,m)=(pow(-1.0,m+1))*adl;
    }
    
    //#% ! Now loop over the other columns backwards
    
    for(int k = m-1;k>=0;k--){
        for(int i = 0;i<=2*m+1;i++){
            double adl, adr, c1l, c1r, c2l, c2r;
            adl = 0.0;
            adr = 0.0;
            for(int j = std::max(0,i-k);j<=std::min(i,m+1);j++){
                if((k-i+j)==0){
                    c2l = 1.0;
                    c2r = 1.0;
                }
                else{
                    c2l = pow(z,k-i+j);
                    c2r = pow(zc,k-i+j);
                }
                if((m+1-j)==0){
                    c1l = 1.0;
                    c1r = 1.0;
                }
                else{
                    c1l = pow(zc,m+1-j);
                    c1r = pow(z,m+1-j);
                }
                adr += bcofs(m+1,j)*bcofs(k,i-j)*c1r*c2r;
                adl += bcofs(m+1,j)*bcofs(k,i-j)*c1l*c2l;
            }
            tmat(i,k+m+1)=adr;
            tmat(i,k)=(pow(-1.0,m+1))*adl;
            double sign;
            sign=1.0;
            for(int j = k+1;j<=m;j++){
                sign=-sign;
                tmat(i,k) -= sign*bcofs(m+1,j-k)*tmat(i,j);
                tmat(i,k+m+1) -= bcofs(m+1,j-k)*tmat(i,j+m+1);
            }
        }
    }
}


//------------------------------------------
void Hermite::interpolate_pd(Darray4 &u, Darray4 &v){
    // This subroutine interpolates from primal to dual

    // LAPACK info
    int M,N,LDA,ONE;
    int Mv,Nv,LDAv;
    char no;
    double ALPHA,BETA;
    M = N = LDA = 2*m+2;
    Mv = Nv = LDAv = 2*m;
    ONE = 1;
    no = 'N';
    ALPHA = 1.0;
    BETA = 0.0;

    // Create pointers to the maps for LAPACK routines
    double *Hmap_u = Hmat_u.c_ptr();
    double *Hmap_v = Hmat_v.c_ptr();

    // Interpolation Arrays
    Darray1 uBottom, uTop;
    uBottom.define(0,2*m+1);
    uBottom.set_value(0.0);
    uTop.define(0,2*m+1);
    uTop.set_value(0.0);

    Darray1 vBottom, vTop;
    vBottom.define(0,2*m-1);
    vBottom.set_value(0.0);
    vTop.define(0,2*m-1);
    vTop.set_value(0.0);

    Darray1 uBottom_int, uTop_int;
    uBottom_int.define(0,2*m+1);
    uBottom_int.set_value(0.0);
    uTop_int.define(0,2*m+1);
    uTop_int.set_value(0.0);

    Darray1 vBottom_int, vTop_int;
    vBottom_int.define(0,2*m-1);
    vBottom_int.set_value(0.0);
    vTop_int.define(0,2*m-1);
    vTop_int.set_value(0.0);

    Darray2 uBottom_2D, uTop_2D;
    uBottom_2D.define(0,2*m+1,0,m);
    uBottom_2D.set_value(0.0);
    uTop_2D.define(0,2*m+1,0,m);
    uTop_2D.set_value(0.0);

    Darray2 vBottom_2D, vTop_2D;
    vBottom_2D.define(0,2*m-1,0,m);
    vBottom_2D.set_value(0.0);
    vTop_2D.define(0,2*m-1,0,m);
    vTop_2D.set_value(0.0);

    // Array for solution
    ud_interp.define(0,2*m+1,0,2*m+1,1,nr,1,ns);
    ud_interp.set_value(0.0);
    vd_interp.define(0,2*m-1,0,2*m-1,1,nr,1,ns);
    vd_interp.set_value(0.0);




    // Interpolate onto dual grid
    for (int i = 0; i <= nr-1; i++){
        for (int j = 0; j <= ns-1; j++){
            for (int ids = 0; ids <= m; ids++){
                for (int idr = 0; idr <= m; idr++){
                    uBottom(idr) = u(idr,ids,i,j);
                    uBottom(m+1+idr) = u(idr,ids,i+1,j);
                    uTop(idr) = u(idr,ids,i,j+1);
                    uTop(m+1+idr) = u(idr,ids,i+1,j+1); 

                }
                // Here we interpolate onto the points (r_i+1/2,s_j) and (r_i+1/2,s_j+1)
                double *uBottom_ptr = uBottom.c_ptr();
                double *uTop_ptr = uTop.c_ptr();
                double *uBottom_int_ptr = uBottom_int.c_ptr();
                double *uTop_int_ptr = uTop_int.c_ptr();
                dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,uBottom_ptr,&ONE,&BETA,uBottom_int_ptr,&ONE);
                dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,uTop_ptr,&ONE,&BETA,uTop_int_ptr,&ONE);
                for(int i=0;i<=2*m+1;i++){
                    uBottom_2D(i,ids) = uBottom_int(i);
                    uTop_2D(i,ids) = uTop_int(i);
                }
            }

            // Now interpolate in y direction
            Darray1 y_interp;
            y_interp.define(0,2*m+1);
            y_interp.set_value(0.0);
            Darray1 interp;
            interp.define(0,2*m+1);
            interp.set_value(0.0);
            double *interp_ptr = interp.c_ptr();
            for(int idr = 0; idr <= 2*m+1; idr++){
                for(int ids = 0; ids <= m; ids++){
                    y_interp(ids) = uBottom_2D(idr,ids);
                    y_interp(m+1+ids) = uTop_2D(idr,ids);
                }
                double *y_interp_ptr = y_interp.c_ptr();
                // Here we interpolate onto the points (r_i+1/2,s_j+1/2)
                dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,y_interp_ptr,&ONE,&BETA,interp_ptr,&ONE);
                for(int ids = 0; ids <= 2*m+1; ids++){
                    ud_interp(idr,ids,i+1,j+1) = interp(ids);
                }
            }

            for (int ids = 0; ids <= m-1; ids++){
                for (int idr = 0; idr <= m-1; idr++){
                    vBottom(idr) = v(idr,ids,i,j);
                    vBottom(m+idr) = v(idr,ids,i+1,j);
                    vTop(idr) = v(idr,ids,i,j+1);
                    vTop(m+idr) = v(idr,ids,i+1,j+1); 

                }
                // Here we interpolate onto the points (r_i+1/2,s_j) and (r_i+1/2,s_j+1)
                double *vBottom_ptr = vBottom.c_ptr();
                double *vTop_ptr = vTop.c_ptr();
                double *vBottom_int_ptr = vBottom_int.c_ptr();
                double *vTop_int_ptr = vTop_int.c_ptr();
                dgemv_(&no,&Mv,&Nv,&ALPHA,Hmap_v,&Nv,vBottom_ptr,&ONE,&BETA,vBottom_int_ptr,&ONE);
                dgemv_(&no,&Mv,&Nv,&ALPHA,Hmap_v,&Nv,vTop_ptr,&ONE,&BETA,vTop_int_ptr,&ONE);
                for(int i=0;i<=2*m-1;i++){
                    vBottom_2D(i,ids) = vBottom_int(i);
                    vTop_2D(i,ids) = vTop_int(i);
                }
            }

            // Now interpolate in y direction
            Darray1 y_interp_v;
            y_interp_v.define(0,2*m-1);
            y_interp_v.set_value(0.0);
            Darray1 interp_v;
            interp_v.define(0,2*m-1);
            interp_v.set_value(0.0);
            double *interp_v_ptr = interp_v.c_ptr();
            for(int idr = 0; idr <= 2*m-1; idr++){
                for(int ids = 0; ids <= m-1; ids++){
                    y_interp_v(ids) = vBottom_2D(idr,ids);
                    y_interp_v(m+ids) = vTop_2D(idr,ids);
                }
                double *y_interp_v_ptr = y_interp_v.c_ptr();
                // Here we interpolate onto the points (r_i+1/2,s_j+1/2)
                dgemv_(&no,&Mv,&Nv,&ALPHA,Hmap_v,&Nv,y_interp_v_ptr,&ONE,&BETA,interp_v_ptr,&ONE);
                for(int ids = 0; ids <= 2*m-1; ids++){
                    vd_interp(idr,ids,i+1,j+1) = interp_v(ids);
                }
            }
        }
    }
}

//--------------------------------------------------
void Hermite::interpolate_dp(Darray4 &ud, Darray4 &vd){
    // This subroutine interpolates from dual to priaml

    // LAPACK info
    int M,N,LDA,ONE;
    int Mv,Nv,LDAv;
    char no;
    double ALPHA,BETA;
    M = N = LDA = 2*m+2;
    Mv = Nv = LDAv = 2*m;
    ONE = 1;
    no = 'N';
    ALPHA = 1.0;
    BETA = 0.0;

    // Create pointers to the maps for LAPACK routines
    double *Hmap_u = Hmat_u.c_ptr();
    double *Hmap_v = Hmat_v.c_ptr();

    // Interpolation Arrays
    Darray1 uBottom, uTop;
    uBottom.define(0,2*m+1);
    uBottom.set_value(0.0);
    uTop.define(0,2*m+1);
    uTop.set_value(0.0);

    Darray1 vBottom, vTop;
    vBottom.define(0,2*m-1);
    vBottom.set_value(0.0);
    vTop.define(0,2*m-1);
    vTop.set_value(0.0);

    Darray1 uBottom_int, uTop_int;
    uBottom_int.define(0,2*m+1);
    uBottom_int.set_value(0.0);
    uTop_int.define(0,2*m+1);
    uTop_int.set_value(0.0);

    Darray1 vBottom_int, vTop_int;
    vBottom_int.define(0,2*m-1);
    vBottom_int.set_value(0.0);
    vTop_int.define(0,2*m-1);
    vTop_int.set_value(0.0);

    Darray2 uBottom_2D, uTop_2D;
    uBottom_2D.define(0,2*m+1,0,m);
    uBottom_2D.set_value(0.0);
    uTop_2D.define(0,2*m+1,0,m);
    uTop_2D.set_value(0.0);

    Darray2 vBottom_2D, vTop_2D;
    vBottom_2D.define(0,2*m-1,0,m);
    vBottom_2D.set_value(0.0);
    vTop_2D.define(0,2*m-1,0,m);
    vTop_2D.set_value(0.0);


    // Interpolate onto dual grid
    for (int i = 1; i <= nr-1; i++){
        for (int j = 1; j <= ns-1; j++){
            for (int ids = 0; ids <= m; ids++){
                for (int idr = 0; idr <= m; idr++){
                    uBottom(idr) = ud(idr,ids,i,j);
                    uBottom(m+1+idr) = ud(idr,ids,i+1,j);
                    uTop(idr) = ud(idr,ids,i,j+1);
                    uTop(m+1+idr) = ud(idr,ids,i+1,j+1); 

                }
                // Here we interpolate onto the points (r_i+1/2,s_j) and (r_i+1/2,s_j+1)
                double *uBottom_ptr = uBottom.c_ptr();
                double *uTop_ptr = uTop.c_ptr();
                double *uBottom_int_ptr = uBottom_int.c_ptr();
                double *uTop_int_ptr = uTop_int.c_ptr();
                dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,uBottom_ptr,&ONE,&BETA,uBottom_int_ptr,&ONE);
                dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,uTop_ptr,&ONE,&BETA,uTop_int_ptr,&ONE);
                for(int i=0;i<=2*m+1;i++){
                    uBottom_2D(i,ids) = uBottom_int(i);
                    uTop_2D(i,ids) = uTop_int(i);
                }
            }

            // Now interpolate in y direction
            Darray1 y_interp;
            y_interp.define(0,2*m+1);
            y_interp.set_value(0.0);
            Darray1 interp;
            interp.define(0,2*m+1);
            interp.set_value(0.0);
            double *interp_ptr = interp.c_ptr();
            for(int idr = 0; idr <= 2*m+1; idr++){
                for(int ids = 0; ids <= m; ids++){
                    y_interp(ids) = uBottom_2D(idr,ids);
                    y_interp(m+1+ids) = uTop_2D(idr,ids);
                }
                double *y_interp_ptr = y_interp.c_ptr();
                // Here we interpolate onto the points (r_i+1/2,s_j+1/2)
                dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,y_interp_ptr,&ONE,&BETA,interp_ptr,&ONE);
                for(int ids = 0; ids <= 2*m+1; ids++){
                    u_interp(idr,ids,i,j) = interp(ids);
                }
            }

            for (int ids = 0; ids <= m-1; ids++){
                for (int idr = 0; idr <= m-1; idr++){
                    vBottom(idr) = vd(idr,ids,i,j);
                    vBottom(m+idr) = vd(idr,ids,i+1,j);
                    vTop(idr) = vd(idr,ids,i,j+1);
                    vTop(m+idr) = vd(idr,ids,i+1,j+1); 

                }
                // Here we interpolate onto the points (r_i+1/2,s_j) and (r_i+1/2,s_j+1)
                double *vBottom_ptr = vBottom.c_ptr();
                double *vTop_ptr = vTop.c_ptr();
                double *vBottom_int_ptr = vBottom_int.c_ptr();
                double *vTop_int_ptr = vTop_int.c_ptr();
                dgemv_(&no,&Mv,&Nv,&ALPHA,Hmap_v,&Nv,vBottom_ptr,&ONE,&BETA,vBottom_int_ptr,&ONE);
                dgemv_(&no,&Mv,&Nv,&ALPHA,Hmap_v,&Nv,vTop_ptr,&ONE,&BETA,vTop_int_ptr,&ONE);
                for(int i=0;i<=2*m-1;i++){
                    vBottom_2D(i,ids) = vBottom_int(i);
                    vTop_2D(i,ids) = vTop_int(i);
                }
            }

            // Now interpolate in y direction
            Darray1 y_interp_v;
            y_interp_v.define(0,2*m-1);
            y_interp_v.set_value(0.0);
            Darray1 interp_v;
            interp_v.define(0,2*m-1);
            interp_v.set_value(0.0);
            double *interp_v_ptr = interp_v.c_ptr();
            for(int idr = 0; idr <= 2*m-1; idr++){
                for(int ids = 0; ids <= m-1; ids++){
                    y_interp_v(ids) = vBottom_2D(idr,ids);
                    y_interp_v(m+ids) = vTop_2D(idr,ids);
                }
                double *y_interp_v_ptr = y_interp_v.c_ptr();
                // Here we interpolate onto the points (r_i+1/2,s_j+1/2)
                dgemv_(&no,&Mv,&Nv,&ALPHA,Hmap_v,&Nv,y_interp_v_ptr,&ONE,&BETA,interp_v_ptr,&ONE);
                for(int ids = 0; ids <= 2*m-1; ids++){
                    v_interp(idr,ids,i,j) = interp_v(ids);
                }
            }
        }
    }
} 

//------------------------------------------
void Hermite::cartesian_recursion(Darray4 &u_int,Darray4 &v_int,Darray4 &u,Darray4 &v,int i,int j,double hr, double hs, double dt){
    int m2p1 = 2*m+1;
    int q = 4*m-1;

    // Create U and V arrays in recursion
    Darray3 U, V;
    U.define(0,2*m+1,0,2*m+1,0,q);
    V.define(0,2*m+1,0,2*m+1,0,q);
    U.set_value(0.0);
    V.set_value(0.0);

    // Set s = 0
    for(int k = 0; k <= 2*m+1; k++){
        for(int l = 0; l <= 2*m+1; l++){
            U(k,l,0) = u_int(k,l,i,j);
        }
    }
    for(int k = 0; k <= 2*m-1; k++){
        for(int l = 0; l<=2*m-1; l++){
            V(k,l,0) = v_int(k,l,i,j);
        }
    }
    // Cartesian Recursion 
    for (int s = 1; s <= q; s++){
        for (int l = 0; l <= 2*m-1; l++){
            for (int k = 0; k <= 2*m-1; k++){
                U(k,l,s) = (dt/double(s))*V(k,l,s-1); 
                V(k,l,s) = dt*(k+2)*(k+1)*U(k+2,l,s-1)/(hr*hr*s) + dt*(l+2)*(l+1)*U(k,l+2,s-1)/(hs*hs*s);
            }
        }
    }

    // Update 
    for(int l = 0; l <= m; l++){
        for(int k = 0; k <=m; k++){
            double Uval = 0.0;
            for(int s =0; s <= q; s++){
                Uval += U(k,l,s)*pow(0.5,s);
            }
            u(k,l,i,j) = Uval;
        }
    }
    for(int l = 0; l <= m-1; l++){
        for(int k = 0; k <= m-1; k++){
            double Vval = 0.0;
            for(int s = 0; s <= q; s++){
                Vval += V(k,l,s)*pow(0.5,s);
            }
            v(k,l,i,j) = Vval;
        }
    }
}

//--------------------------------------------------


void Hermite::recursion_curvilinear(Darray4 &u_int,Darray4 &v_int,int i,int j,
    Darray4 &A01, Darray4 &A10,Darray4 & A11,Darray4 &A12,Darray4 & A22,Darray4 &u,Darray4 &v,
    double hr,double hs,double dt){
    int m2p1 = 2*m+1;
    int q = 4*m-1;
    Darray2 ur,us,urs,urr,uss;
    ur.define(0,m2p1,0,m2p1);
    us.define(0,m2p1,0,m2p1);
    urs.define(0,m2p1,0,m2p1);
    urr.define(0,m2p1,0,m2p1);
    uss.define(0,m2p1,0,m2p1);
    ur.set_value(0.0);
    us.set_value(0.0);
    urs.set_value(0.0);
    urr.set_value(0.0);
    uss.set_value(0.0);
    double *ur_ptr = ur.c_ptr();
    double *us_ptr = us.c_ptr();
    double *urs_ptr = urs.c_ptr();
    double *urr_ptr = urr.c_ptr();
    double *uss_ptr = uss.c_ptr();
    // Define arrays needed for computation
    Darray2 lap;
    lap.define(0,m2p1,0,m2p1);
    lap.set_value(0.0);
    double *lap_ptr = lap.c_ptr();
    // Define variable coefficients for recursion
    Darray2 a01_recursion, a10_recursion, a11_recursion, a12_recursion, a22_recursion;
    a01_recursion.define(0,m2p1,0,m2p1);
    a10_recursion.define(0,m2p1,0,m2p1);
    a11_recursion.define(0,m2p1,0,m2p1);
    a12_recursion.define(0,m2p1,0,m2p1);
    a22_recursion.define(0,m2p1,0,m2p1);
    a01_recursion.set_value(0.0);
    a10_recursion.set_value(0.0);
    a11_recursion.set_value(0.0);
    a12_recursion.set_value(0.0);
    a22_recursion.set_value(0.0);

    // Pointers 
    double* a01_recursion_ptr = a01_recursion.c_ptr();
    double* a10_recursion_ptr = a10_recursion.c_ptr();
    double* a11_recursion_ptr = a11_recursion.c_ptr();
    double* a12_recursion_ptr = a12_recursion.c_ptr();
    double* a22_recursion_ptr = a22_recursion.c_ptr();

    
    // Set variable coefficient values
    for(int idr = 0; idr <= m2p1; idr++){
        for(int ids = 0; ids <= m2p1; ids++){
            a01_recursion(idr,ids) = A01(idr,ids,i,j);
            a10_recursion(idr,ids) = A10(idr,ids,i,j);
            a11_recursion(idr,ids) = A11(idr,ids,i,j);
            a12_recursion(idr,ids) = A12(idr,ids,i,j);
            a22_recursion(idr,ids) = A22(idr,ids,i,j);
        }
    }

    Darray2 a11_urr, a12_urs, a22_uss, a10_ur, a01_us;
    a11_urr.define(0,m2p1,0,m2p1);
    a12_urs.define(0,m2p1,0,m2p1);
    a22_uss.define(0,m2p1,0,m2p1);
    a10_ur.define(0,m2p1,0,m2p1);
    a01_us.define(0,m2p1,0,m2p1);
    a11_urr.set_value(0.0);
    a12_urs.set_value(0.0);
    a22_uss.set_value(0.0);
    a10_ur.set_value(0.0);
    a01_us.set_value(0.0);

    double* a11_urr_ptr = a11_urr.c_ptr();
    double* a12_urs_ptr = a12_urs.c_ptr();
    double* a22_uss_ptr = a22_uss.c_ptr();
    double* a10_ur_ptr = a10_ur.c_ptr();
    double* a01_us_ptr = a01_us.c_ptr();

    // Create U and V arrays in recursion
    Darray3 U, V;
    U.define(0,m2p1,0,m2p1,0,q);
    V.define(0,m2p1,0,m2p1,0,q);
    U.set_value(0.0);
    V.set_value(0.0);
    for(int k = 0; k <= m2p1; k++){
        for(int l = 0; l <= m2p1; l++){
            U(k,l,0) = u_int(k,l,i,j);
        }
    }
    for(int k = 0; k <= 2*m-1; k++){
        for(int l = 0; l<=2*m-1; l++){
            V(k,l,0) = v_int(k,l,i,j);
        }
    }
    // Time Derivatives
    for(int s = 1; s <= q; s++){
        for(int l = 0; l <= 2*m-1; l++){
            for(int k = 0; k <= 2*m-1; k++){
                U(k,l,s) = dt*V(k,l,s-1)/double(s);
            }
        }
        for(int l = 0; l <= 2*m+1; l++){
            for(int k = 0; k <= 2*m; k++){
                ur(k,l) = double(k+1)*U(k+1,l,s-1)/hr;
            }
        }
        for(int l = 0; l <= 2*m; l++){
            for(int k = 0; k <= 2*m+1; k++){
                us(k,l) = double(l+1)*U(k,l+1,s-1)/hs;
            }
        }
        for(int l = 0; l <= 2*m+1; l++){
            for(int k = 0; k <= 2*m-1; k++){
                urr(k,l) = double(k+1)*double(k+2)*U(k+2,l,s-1)/(hr*hr);
            }
        }
        for(int l = 0; l <= 2*m-1; l++){
            for(int k = 0; k <= 2*m+1; k++){
                uss(k,l) = double(l+1)*double(l+2)*U(k,l+2,s-1)/(hs*hs);
            }
        }
        for(int l = 0; l <= 2*m; l++){
            for(int k = 0; k <= 2*m; k++){
                urs(k,l) = double(l+1)*double(k+1)*U(k+1,l+1,s-1)/(hr*hs);
            }
        }

        tpolymul2_(a11_recursion_ptr,urr_ptr,a11_urr_ptr,&m2p1);
        tpolymul2_(a12_recursion_ptr,urs_ptr,a12_urs_ptr,&m2p1);
        tpolymul2_(a22_recursion_ptr,uss_ptr,a22_uss_ptr,&m2p1);
        tpolymul2_(a10_recursion_ptr,ur_ptr,a10_ur_ptr,&m2p1);
        tpolymul2_(a01_recursion_ptr,us_ptr,a01_us_ptr,&m2p1);
        for(int idr = 0; idr <= m2p1; idr++){
            for(int ids = 0; ids <= m2p1; ids++){
                lap(idr,ids) = a11_urr(idr,ids) + a12_urs(idr,ids) + a22_uss(idr,ids)
                + a10_ur(idr,ids) + a01_us(idr,ids);
            }
        }
        for(int idr = 0; idr <= m2p1-2; idr++){
            for(int ids = 0; ids <= m2p1-2; ids++){
                V(idr,ids,s) = dt*lap(idr,ids)/double(s);
            }
        }
        lap.set_value(0.0);
    }
    for(int idr = 0; idr <= m; idr++){
        for(int ids = 0; ids <= m; ids++){
            double Uval = 0.0;
            for(int s = 0; s <= q; s++){
                Uval += U(idr,ids,s)*pow(0.5,s);
            }
            u(idr,ids,i,j) = Uval;
        }
    }
    for(int idr = 0; idr <= m-1; idr++){
        for(int ids = 0; ids <= m-1; ids++){
            double Vval = 0.0;
            for(int s = 0; s <= q; s++){
                Vval += V(idr,ids,s)*pow(0.5,s);
            }
            v(idr,ids,i,j) = Vval;
        }
    }
}

void Hermite::outPutMat(Darray2 &M,int i_start,int i_end,int j_start,int j_end,char* fName){
    FILE *extFile = fopen(fName, "w");
    for (int i=i_start; i<=i_end; i++){
        for (int j=j_start; j <= j_end; j++){
            fprintf(extFile, " %18.10e ", M(i,j));
        }
        fprintf(extFile,"\n");
    }   
    fclose(extFile);
}
void Hermite::outPutVec(Darray1 &M,int i_start,int i_end,char* fName){
    FILE *extFile = fopen(fName, "w");
    for (int i=i_start; i<=i_end; i++){
            fprintf(extFile, "%18.10e \n ", M(i));
        }
    fclose(extFile);
}

void columnScale(Darray2 &A,Darray1 &b, int m){
    int matSize = (2*m+2)*(2*m+2)-1;
    Darray1 cSum;
    cSum.define(0,matSize);
    cSum.set_value(0.0);

    double colVal = 0.0;
    for(int i = 0; i <= matSize; i++){
        for(int j = 0; j <= matSize; j++){
            colVal += abs(A(i,j));
        }
        cSum(i) = colVal;
        colVal = 0.0;
    }

    for(int j = 0; j <= matSize; j++){
        for(int i = 0; i <= matSize; i++){
            A(i,j) = A(i,j)/cSum(i);
        }
    }

    for(int j = 0; j <= matSize; j++){
        b(j) = b(j)/cSum(j);
        // cout << j << " cSum " << cSum(j) << endl;
    }
}
void Hermite::computeSingular(int m,Darray2 Mat_u,int ns){

        // SVD LAPACK STUFF
        int Mboundary = (2*m+2)*(2*m+2);
        int INFO;
        char JOBU = 'A';
        char JOBVT = 'A';
        int LDA = (2*m+2)*(2*m+2);
        int LDVT = (2*m+2)*(2*m+2);
        int LDU = (2*m+2)*(2*m+2);
        int Lwork = 5*((2*m+2)*(2*m+2));
        Darray1 S;
        S.define(0,(2*m+2)*(2*m+2)-1);
        S.set_value(0.0);
        double * S_ptr = S.c_ptr();
        Darray2 U, VT,work;
        U.define(0,(2*m+2)*(2*m+2)-1,0,(2*m+2)*(2*m+2)-1);
        VT.define(0,(2*m+2)*(2*m+2)-1,0,(2*m+2)*(2*m+2)-1);
        work.define(0,(2*m+2)*(2*m+2)-1,0,(2*m+2)*(2*m+2)-1);
        U.set_value(0.0);
        VT.set_value(0.0);
        work.set_value(0.0);
        double * U_ptr = U.c_ptr();
        double * VT_ptr = VT.c_ptr();
        double * work_ptr = work.c_ptr();
        // Define arrays
        // Darray2 a10, a01, a11, a12, a22;
        // a10.define(0,2*m+1+2,0,2*m+1+2);
        // a01.define(0,2*m+1+2,0,2*m+1+2);
        // a11.define(0,2*m+1+2,0,2*m+1+2);
        // a12.define(0,2*m+1+2,0,2*m+1+2);
        // a22.define(0,2*m+1+2,0,2*m+1+2);
        // a10.set_value(0.0);
        // a01.set_value(0.0);
        // a11.set_value(0.0);
        // a12.set_value(0.0);
        // a22.set_value(0.0);
        // double* a10_ptr = a10.c_ptr();
        // double* a01_ptr = a01.c_ptr();
        // double* a11_ptr = a11.c_ptr();
        // double* a12_ptr = a12.c_ptr();
        // double* a22_ptr = a22.c_ptr();
        // Darray2 Mat_u;
        // Mat_u.define(0,(2*m+2)*(2*m+2)-1,0,(2*m+2)*(2*m+2)-1);
        // Mat_u.set_value(0.0);
        double* Mat_u_ptr = Mat_u.c_ptr();
        // Create matrix
        // double rA,rB,sA,sB;
        // rA = 0.5;
        // rB = 0.5;
        // sA = 0.5;
        // sB = -0.5;
        // int m2p1 = 2*m+1;
        // // Set variable Coefficient values
        // for(int idr = 0; idr <= m2p1; idr++){
        //     for(int ids = 0; ids <= m2p1; ids++){
        //         a01(idr,ids) = A01(idr,ids,0,1)/(pow(hr,idr)*pow(hs,ids));
        //         a10(idr,ids) = A10(idr,ids,0,1)/(pow(hr,idr)*pow(hs,ids));
        //         a11(idr,ids) = A11(idr,ids,0,1)/(pow(hr,idr)*pow(hs,ids));
        //         a12(idr,ids) = A12(idr,ids,0,1)/(pow(hr,idr)*pow(hs,ids));
        //         a22(idr,ids) = A22(idr,ids,0,1)/(pow(hr,idr)*pow(hs,ids));
        //     }
        // }

    // bcMat(rA,rB,sA,sB,hr,hs,a01_ptr,a10_ptr,a11_ptr,a12_ptr,a22_ptr,Mat_u_ptr,m);
    dgesvd_(& JOBU,& JOBVT,&Mboundary,&Mboundary,Mat_u_ptr,&Mboundary,S_ptr,U_ptr,&Mboundary,VT_ptr,&Mboundary,work_ptr,&Lwork,&INFO);
    cout << "INFO = " << INFO << endl;
    char fName[100];
    sprintf(fName, "Singular%0.5d.ext",ns);
    outPutVec(S,0,(2*m+2)*(2*m+2)-1,fName);
}
void Hermite::boundaryConditions(int m,double hr,double hs,Darray4 &u,Darray4 &v, 
    Darray4 &A01,Darray4 & A10,Darray4 & A11,Darray4 &A12,Darray4 & A22){
    int m2p1 = 2*m+1;
    int m_minus = m - 1;
    // LAPACK info for u
    int Mboundary,Nboundary,LDAboundary;
    Mboundary = Nboundary = LDAboundary = (2*m+2)*(2*m+2);
    int IPIV[(2*m+2)*(2*m+2)];
    // LAPACK info for v
    int Mboundary_v,Nboundary_v,LDAboundary_v;
    Mboundary_v = Nboundary_v = LDAboundary_v = (2*m)*(2*m);
    int IPIV_v[(2*m)*(2*m)];
    // LAPACK info for both
    int INFO;
    int ONE = 1;
    char no = 'N';
    char fName[100];

    // Define arrays
    Darray2 a10, a01, a11, a12, a22;
    a10.define(0,2*m+1+2,0,2*m+1+2);
    a01.define(0,2*m+1+2,0,2*m+1+2);
    a11.define(0,2*m+1+2,0,2*m+1+2);
    a12.define(0,2*m+1+2,0,2*m+1+2);
    a22.define(0,2*m+1+2,0,2*m+1+2);
    a10.set_value(0.0);
    a01.set_value(0.0);
    a11.set_value(0.0);
    a12.set_value(0.0);
    a22.set_value(0.0);
    double* a10_ptr = a10.c_ptr();
    double* a01_ptr = a01.c_ptr();
    double* a11_ptr = a11.c_ptr();
    double* a12_ptr = a12.c_ptr();
    double* a22_ptr = a22.c_ptr();
    Darray2 Mat_u;
    Mat_u.define(0,(2*m+2)*(2*m+2)-1,0,(2*m+2)*(2*m+2)-1);
    Mat_u.set_value(0.0);
    double* Mat_u_ptr = Mat_u.c_ptr();
    Darray2 Mat_v;
    Mat_v.define(0,(2*m)*(2*m)-1,0,(2*m)*(2*m)-1);
    Mat_v.set_value(0.0);
    double* Mat_v_ptr = Mat_v.c_ptr();
    // Declare parameters
    double rA,rB,sA,sB;
    Darray1 rhs_u;
    rhs_u.define(0,(2*m+2)*(2*m+2)-1);
    double* rhs_u_ptr = rhs_u.c_ptr();
    Darray1 rhs_v;
    rhs_v.define(0,(2*m)*(2*m)-1);
    double* rhs_v_ptr = rhs_v.c_ptr();
    int offset = 0;
    int COUNT = 0;

    Darray1 rhs_u_perm;
    rhs_u_perm.define(0,(2*m+2)*(2*m+2)-1);
    rhs_u_perm.set_value(0.0);

    Darray1 rhs_v_perm;
    rhs_v_perm.define(0,(2*m)*(2*m)-1);
    rhs_v_perm.set_value(0.0);

    //  Equilibrium
    Darray1 A;
    A.define(0,(2*m+2)*(2*m+2)*(2*m+2)*(2*m+2)-1);
    A.set_value(0.0);
    double *A_ptr = A.c_ptr();

    int N, NE,LENC_number;
    N = (2*m+2)*(2*m+2);
    NE = (2*m+2)*(2*m+2)*(2*m+2)*(2*m+2);
    LENC_number = (2*m+2)*(2*m+2);
    Darray2 A_mat;
    A_mat.define(0,(2*m+2)*(2*m+2)-1,0,(2*m+2)*(2*m+2)-1);
    A_mat.set_value(0.0);
    double *A_mat_ptr = A_mat.c_ptr();

    Darray2 B_mat;
    B_mat.define(0,(2*m+2)*(2*m+2)-1,0,(2*m+2)*(2*m+2)-1);
    B_mat.set_value(0.0);
    double *B_mat_ptr = B_mat.c_ptr();

    Darray1 ROW;
    ROW.define(0,(2*m+2)*(2*m+2)-1);
    ROW.set_value(0.0);
    double * ROW_ptr = ROW.c_ptr();

    Darray1 COL;
    COL.define(0,(2*m+2)*(2*m+2)-1);
    COL.set_value(0.0);
    double * COL_ptr = COL.c_ptr();

    Darray1 PERM;
    PERM.define(0,(2*m+2)*(2*m+2)-1);
    PERM.set_value(0.0);
    double * PERM_ptr = PERM.c_ptr();

    // For v
    Darray1 Av;
    Av.define(0,(2*m)*(2*m)*(2*m)*(2*m)-1);
    Av.set_value(0.0);
    double *Av_ptr = Av.c_ptr();

    int Nv, NEv,LENC_numberv;
    Nv = (2*m)*(2*m);
    NEv = (2*m)*(2*m)*(2*m)*(2*m);
    LENC_numberv = (2*m)*(2*m);
    Darray2 A_matv;
    A_matv.define(0,(2*m)*(2*m)-1,0,(2*m)*(2*m)-1);
    A_matv.set_value(0.0);
    double *A_matv_ptr = A_matv.c_ptr();

    Darray2 B_matv;
    B_matv.define(0,(2*m)*(2*m)-1,0,(2*m)*(2*m)-1);
    B_matv.set_value(0.0);
    double *B_matv_ptr = B_matv.c_ptr();

    Darray1 ROWv;
    ROWv.define(0,(2*m)*(2*m)-1);
    ROWv.set_value(0.0);
    double * ROWv_ptr = ROWv.c_ptr();

    Darray1 COLv;
    COLv.define(0,(2*m)*(2*m)-1);
    COLv.set_value(0.0);
    double * COLv_ptr = COLv.c_ptr();

    Darray1 PERMv;
    PERMv.define(0,(2*m)*(2*m)-1);
    PERMv.set_value(0.0);
    double * PERMv_ptr = PERMv.c_ptr();

    for(int j = 1; j <= ns - 1; j++){
        // Left side
        rA = 0.5;
        rB = 0.5;
        sA = 0.5;
        sB = -0.5;
        // Set variable Coefficient values
        for(int idr = 0; idr <= m2p1; idr++){
            for(int ids = 0; ids <= m2p1; ids++){
                a01(idr,ids) = A01(idr,ids,0,j)/(pow(hr,idr)*pow(hs,ids));
                a10(idr,ids) = A10(idr,ids,0,j)/(pow(hr,idr)*pow(hs,ids));
                a11(idr,ids) = A11(idr,ids,0,j)/(pow(hr,idr)*pow(hs,ids));
                a12(idr,ids) = A12(idr,ids,0,j)/(pow(hr,idr)*pow(hs,ids));
                a22(idr,ids) = A22(idr,ids,0,j)/(pow(hr,idr)*pow(hs,ids));
            }
        }
        // Boundary Conditions for u
        rhs_u.set_value(0.0);
        offset = 0;
        for(int ids = 0;ids<=m;ids++){
            for(int idr = 0;idr<=m;idr++){
                rhs_u((2*m+2)*(2*m+2)/2+offset) = tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,1,j+1)/(pow(hr,idr)*pow(hs,ids));
                offset++;
                rhs_u((2*m+2)*(2*m+2)/2+offset) =  tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,1,j)/(pow(hr,idr)*pow(hs,ids));
                offset++;
            }
        }
        bcMat(rA,rB,sA,sB,hr,hs,a01_ptr,a10_ptr,a11_ptr,a12_ptr,a22_ptr,Mat_u_ptr,m);
        // Equilibriate
        COUNT = 0;
        for(int j = 0; j <= (2*m+2)*(2*m+2)-1;j++){
            for(int i = 0; i <= (2*m+2)*(2*m+2)-1;i++){
                A(COUNT) = Mat_u(i,j);      
                COUNT += 1;
            }
        }
        eq_(A_ptr,&N,&NE,&LENC_number,ROW_ptr,COL_ptr,PERM_ptr);
        COUNT = 0;
        for(int j = 0; j <= (2*m+2)*(2*m+2)-1;j++){
            for(int i = 0; i <= (2*m+2)*(2*m+2)-1;i++){
                A_mat(i,j) = Mat_u(i,j)*exp(ROW(i)+COL(j));      
            }
        }
        for(int j = 0; j <= (2*m+2)*(2*m+2)-1; j++){
            for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
                B_mat(i,j) = A_mat(i,PERM(j)-1);
            }
        }
        sprintf(fName, "M%0.5d.ext",j);
        outPutMat(Mat_u,0,(2*m+2)*(2*m+2)-1,0,(2*m+2)*(2*m+2)-1,fName);
        sprintf(fName, "B%0.5d.ext",j);
        outPutMat(B_mat,0,(2*m+2)*(2*m+2)-1,0,(2*m+2)*(2*m+2)-1,fName);
        dgetrf_(&Mboundary,&Nboundary,B_mat_ptr,&LDAboundary,IPIV,&INFO); 
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
            rhs_u(i) *= exp(ROW(i));
        }
        dgetrs_(&no,&Nboundary,&ONE,B_mat_ptr,&Nboundary,IPIV,rhs_u_ptr,&Nboundary,&INFO); // Solve system
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
            rhs_u_perm(PERM(i)-1) = rhs_u(i);
        }
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
            rhs_u(i) = rhs_u_perm(i)*exp(COL(i));
        }
        // End Equilibriate

        COUNT = 0;
        for(int ids = 0;ids<=2*m+1;ids++){
            for(int idr = 0;idr<=2*m+1;idr++){
                u_interp(idr,ids,0,j) = rhs_u(COUNT)*pow(hr,idr)*pow(hs,ids);
                COUNT++;
            }
        }

        // Boundary Conditions for v
        rhs_v.set_value(0.0);
        offset = 0;
        for(int ids = 0; ids <= m-1; ids++){
            for(int idr = 0; idr <= m-1; idr++){
                rhs_v((2*m)*(2*m)/2+offset) = tgamma(idr+1)*tgamma(ids+1)*v(idr,ids,1,j+1)/(pow(hr,idr)*pow(hs,ids));
                offset++;
                rhs_v((2*m)*(2*m)/2+offset) =  tgamma(idr+1)*tgamma(ids+1)*v(idr,ids,1,j)/(pow(hr,idr)*pow(hs,ids));
                offset++;
            }
        }
        bcMat(rA,rB,sA,sB,hr,hs,a01_ptr,a10_ptr,a11_ptr,a12_ptr,a22_ptr,Mat_v_ptr,m_minus);
        // Equilibriate
        COUNT = 0;
        for(int j = 0; j <= (2*m)*(2*m)-1;j++){
            for(int i = 0; i <= (2*m)*(2*m)-1;i++){
                Av(COUNT) = Mat_v(i,j);      
                COUNT += 1;
            }
        }
        eq_(Av_ptr,&Nv,&NEv,&LENC_numberv,ROWv_ptr,COLv_ptr,PERMv_ptr);
        COUNT = 0;
        for(int j = 0; j <= (2*m)*(2*m)-1;j++){
            for(int i = 0; i <= (2*m)*(2*m)-1;i++){
                A_matv(i,j) = Mat_v(i,j)*exp(ROWv(i)+COLv(j));      
            }
        }
        for(int j = 0; j <= (2*m)*(2*m)-1; j++){
            for(int i = 0; i <= (2*m)*(2*m)-1; i++){
                B_matv(i,j) = A_matv(i,PERMv(j)-1);
            }
        }
        dgetrf_(&Mboundary_v,&Nboundary_v,B_matv_ptr,&LDAboundary_v,IPIV_v,&INFO); 
        for(int i = 0; i <= (2*m)*(2*m)-1; i++){
            rhs_v(i) *= exp(ROWv(i));
        }
        dgetrs_(&no,&Nboundary_v,&ONE,B_matv_ptr,&Nboundary_v,IPIV_v,rhs_v_ptr,&Nboundary_v,&INFO); // Solve system
        for(int i = 0; i <= (2*m)*(2*m)-1; i++){
            rhs_v_perm(PERMv(i)-1) = rhs_v(i);
        }
        for(int i = 0; i <= (2*m)*(2*m)-1; i++){
            rhs_v(i) = rhs_v_perm(i)*exp(COLv(i));
        }
        // End Equilibriate
        COUNT = 0;
        for(int ids = 0;ids<=2*m-1;ids++){
            for(int idr = 0;idr<=2*m-1;idr++){
                v_interp(idr,ids,0,j) = rhs_v(COUNT)*pow(hr,idr)*pow(hs,ids);
                COUNT++;
            }
        }

        // Right side
        rA = -0.5;
        rB = -0.5;
        sA = 0.5;
        sB = -0.5;
        // Set variable Coefficient values
        for(int idr = 0; idr <= m2p1; idr++){
            for(int ids = 0; ids <= m2p1; ids++){
                a01(idr,ids) = A01(idr,ids,nr,j)/(pow(hr,idr)*pow(hs,ids));
                a10(idr,ids) = A10(idr,ids,nr,j)/(pow(hr,idr)*pow(hs,ids));
                a11(idr,ids) = A11(idr,ids,nr,j)/(pow(hr,idr)*pow(hs,ids));
                a12(idr,ids) = A12(idr,ids,nr,j)/(pow(hr,idr)*pow(hs,ids));
                a22(idr,ids) = A22(idr,ids,nr,j)/(pow(hr,idr)*pow(hs,ids));
            }
        }
        // Boundary Conditions for u
        rhs_u.set_value(0.0);
        offset = 0;
        for(int ids = 0;ids<=m;ids++){
            for(int idr = 0;idr<=m;idr++){
                rhs_u((2*m+2)*(2*m+2)/2+offset) = tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,nr,j+1)/(pow(hr,idr)*pow(hs,ids));
                offset++;
                rhs_u((2*m+2)*(2*m+2)/2+offset) =  tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,nr,j)/(pow(hr,idr)*pow(hs,ids));
                offset++;
            }
        }
        bcMat(rA,rB,sA,sB,hr,hs,a01_ptr,a10_ptr,a11_ptr,a12_ptr,a22_ptr,Mat_u_ptr,m);
        // Equilibriate
        COUNT = 0;
        for(int j = 0; j <= (2*m+2)*(2*m+2)-1;j++){
            for(int i = 0; i <= (2*m+2)*(2*m+2)-1;i++){
                A(COUNT) = Mat_u(i,j);      
                COUNT += 1;
            }
        }
        eq_(A_ptr,&N,&NE,&LENC_number,ROW_ptr,COL_ptr,PERM_ptr);
        COUNT = 0;
        for(int j = 0; j <= (2*m+2)*(2*m+2)-1;j++){
            for(int i = 0; i <= (2*m+2)*(2*m+2)-1;i++){
                A_mat(i,j) = Mat_u(i,j)*exp(ROW(i)+COL(j));      
            }
        }
        for(int j = 0; j <= (2*m+2)*(2*m+2)-1; j++){
            for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
                B_mat(i,j) = A_mat(i,PERM(j)-1);
            }
        }
        dgetrf_(&Mboundary,&Nboundary,B_mat_ptr,&LDAboundary,IPIV,&INFO); 
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
            rhs_u(i) *= exp(ROW(i));
        }
        dgetrs_(&no,&Nboundary,&ONE,B_mat_ptr,&Nboundary,IPIV,rhs_u_ptr,&Nboundary,&INFO); // Solve system
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
            rhs_u_perm(PERM(i)-1) = rhs_u(i);
        }
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
            rhs_u(i) = rhs_u_perm(i)*exp(COL(i));
        }
        // End Equilibriate
        COUNT = 0;
        for(int ids = 0;ids<=2*m+1;ids++){
            for(int idr = 0;idr<=2*m+1;idr++){
                u_interp(idr,ids,nr,j) = rhs_u(COUNT)*pow(hr,idr)*pow(hs,ids);
                COUNT++;
            }
        }

        // Boundary Conditions for v
        rhs_v.set_value(0.0);
        offset = 0;
        for(int ids = 0; ids <= m-1; ids++){
            for(int idr = 0; idr <= m-1; idr++){
                rhs_v((2*m)*(2*m)/2+offset) = tgamma(idr+1)*tgamma(ids+1)*v(idr,ids,nr,j+1)/(pow(hr,idr)*pow(hs,ids));
                offset++;
                rhs_v((2*m)*(2*m)/2+offset) =  tgamma(idr+1)*tgamma(ids+1)*v(idr,ids,nr,j)/(pow(hr,idr)*pow(hs,ids));
                offset++;
            }
        }
        bcMat(rA,rB,sA,sB,hr,hs,a01_ptr,a10_ptr,a11_ptr,a12_ptr,a22_ptr,Mat_v_ptr,m_minus);
        // Equilibriate
        COUNT = 0;
        for(int j = 0; j <= (2*m)*(2*m)-1;j++){
            for(int i = 0; i <= (2*m)*(2*m)-1;i++){
                Av(COUNT) = Mat_v(i,j);      
                COUNT += 1;
            }
        }
        eq_(Av_ptr,&Nv,&NEv,&LENC_numberv,ROWv_ptr,COLv_ptr,PERMv_ptr);
        COUNT = 0;
        for(int j = 0; j <= (2*m)*(2*m)-1;j++){
            for(int i = 0; i <= (2*m)*(2*m)-1;i++){
                A_matv(i,j) = Mat_v(i,j)*exp(ROWv(i)+COLv(j));      
            }
        }
        for(int j = 0; j <= (2*m)*(2*m)-1; j++){
            for(int i = 0; i <= (2*m)*(2*m)-1; i++){
                B_matv(i,j) = A_matv(i,PERMv(j)-1);
            }
        }
        dgetrf_(&Mboundary_v,&Nboundary_v,B_matv_ptr,&LDAboundary_v,IPIV_v,&INFO); 
        for(int i = 0; i <= (2*m)*(2*m)-1; i++){
            rhs_v(i) *= exp(ROWv(i));
        }
        dgetrs_(&no,&Nboundary_v,&ONE,B_matv_ptr,&Nboundary_v,IPIV_v,rhs_v_ptr,&Nboundary_v,&INFO); // Solve system
        for(int i = 0; i <= (2*m)*(2*m)-1; i++){
            rhs_v_perm(PERMv(i)-1) = rhs_v(i);
        }
        for(int i = 0; i <= (2*m)*(2*m)-1; i++){
            rhs_v(i) = rhs_v_perm(i)*exp(COLv(i));
        }
        // End Equilibriate
        COUNT = 0;
        for(int ids = 0;ids<=2*m-1;ids++){
            for(int idr = 0;idr<=2*m-1;idr++){
                v_interp(idr,ids,nr,j) = rhs_v(COUNT)*pow(hr,idr)*pow(hs,ids);
                COUNT++;
            }
        }
    }

    // Do corners with periodicity + Boundary conditions
    // Left side
    rA = 0.5;
    rB = 0.5;
    sA = 0.5;
    sB = -0.5;
    // Set variable Coefficient values
    for(int idr = 0; idr <= m2p1; idr++){
        for(int ids = 0; ids <= m2p1; ids++){
            a01(idr,ids) = A01(idr,ids,0,0)/(pow(hr,idr)*pow(hs,ids));
            a10(idr,ids) = A10(idr,ids,0,0)/(pow(hr,idr)*pow(hs,ids));
            a11(idr,ids) = A11(idr,ids,0,0)/(pow(hr,idr)*pow(hs,ids));
            a12(idr,ids) = A12(idr,ids,0,0)/(pow(hr,idr)*pow(hs,ids));
            a22(idr,ids) = A22(idr,ids,0,0)/(pow(hr,idr)*pow(hs,ids));
        }
    }
    // Boundary Conditions for u
    rhs_u.set_value(0.0);
    offset = 0;
    for(int ids = 0;ids<=m;ids++){
        for(int idr = 0;idr<=m;idr++){
            rhs_u((2*m+2)*(2*m+2)/2+offset) = tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,1,1)/(pow(hr,idr)*pow(hs,ids));
            offset++;
            rhs_u((2*m+2)*(2*m+2)/2+offset) =  tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,1,ns)/(pow(hr,idr)*pow(hs,ids));
            offset++;
        }
    }
    bcMat(rA,rB,sA,sB,hr,hs,a01_ptr,a10_ptr,a11_ptr,a12_ptr,a22_ptr,Mat_u_ptr,m);
    columnScale(Mat_u,rhs_u,m);
    // Equilibriate
    COUNT = 0;
    for(int j = 0; j <= (2*m+2)*(2*m+2)-1;j++){
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1;i++){
            A(COUNT) = Mat_u(i,j);      
            COUNT += 1;
        }
    }
    eq_(A_ptr,&N,&NE,&LENC_number,ROW_ptr,COL_ptr,PERM_ptr);
    COUNT = 0;
    for(int j = 0; j <= (2*m+2)*(2*m+2)-1;j++){
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1;i++){
            A_mat(i,j) = Mat_u(i,j)*exp(ROW(i)+COL(j));      
        }
    }
    for(int j = 0; j <= (2*m+2)*(2*m+2)-1; j++){
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
            B_mat(i,j) = A_mat(i,PERM(j)-1);
        }
    }
    dgetrf_(&Mboundary,&Nboundary,B_mat_ptr,&LDAboundary,IPIV,&INFO); 
    for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
        rhs_u(i) *= exp(ROW(i));
    }
    dgetrs_(&no,&Nboundary,&ONE,B_mat_ptr,&Nboundary,IPIV,rhs_u_ptr,&Nboundary,&INFO); // Solve system
    for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
        rhs_u_perm(PERM(i)-1) = rhs_u(i);
    }
    for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
        rhs_u(i) = rhs_u_perm(i)*exp(COL(i));
    }
    // End equilibriate
    COUNT = 0;
    for(int ids = 0;ids<=2*m+1;ids++){
        for(int idr = 0;idr<=2*m+1;idr++){
            u_interp(idr,ids,0,0) = rhs_u(COUNT)*pow(hr,idr)*pow(hs,ids);
            u_interp(idr,ids,0,ns) = rhs_u(COUNT)*pow(hr,idr)*pow(hs,ids); // Use periodicity
            COUNT++;
        }
    }

    // Boundary Conditions for v
    rhs_v.set_value(0.0);
    offset = 0;
    for(int ids = 0; ids <= m-1; ids++){
        for(int idr = 0; idr <= m-1; idr++){
            rhs_v((2*m)*(2*m)/2+offset) = tgamma(idr+1)*tgamma(ids+1)*v(idr,ids,1,1)/(pow(hr,idr)*pow(hs,ids));
            offset++;
            rhs_v((2*m)*(2*m)/2+offset) =  tgamma(idr+1)*tgamma(ids+1)*v(idr,ids,1,ns)/(pow(hr,idr)*pow(hs,ids));
            offset++;
        }
    }
    bcMat(rA,rB,sA,sB,hr,hs,a01_ptr,a10_ptr,a11_ptr,a12_ptr,a22_ptr,Mat_v_ptr,m_minus);
    // Equilibriate
    COUNT = 0;
    for(int j = 0; j <= (2*m)*(2*m)-1;j++){
        for(int i = 0; i <= (2*m)*(2*m)-1;i++){
            Av(COUNT) = Mat_v(i,j);      
            COUNT += 1;
        }
    }
    eq_(Av_ptr,&Nv,&NEv,&LENC_numberv,ROWv_ptr,COLv_ptr,PERMv_ptr);
    COUNT = 0;
    for(int j = 0; j <= (2*m)*(2*m)-1;j++){
        for(int i = 0; i <= (2*m)*(2*m)-1;i++){
            A_matv(i,j) = Mat_v(i,j)*exp(ROWv(i)+COLv(j));      
        }
    }
    for(int j = 0; j <= (2*m)*(2*m)-1; j++){
        for(int i = 0; i <= (2*m)*(2*m)-1; i++){
            B_matv(i,j) = A_matv(i,PERMv(j)-1);
        }
    }
    dgetrf_(&Mboundary_v,&Nboundary_v,B_matv_ptr,&LDAboundary_v,IPIV_v,&INFO); 
    for(int i = 0; i <= (2*m)*(2*m)-1; i++){
        rhs_v(i) *= exp(ROWv(i));
    }
    dgetrs_(&no,&Nboundary_v,&ONE,B_matv_ptr,&Nboundary_v,IPIV_v,rhs_v_ptr,&Nboundary_v,&INFO); // Solve system
    for(int i = 0; i <= (2*m)*(2*m)-1; i++){
        rhs_v_perm(PERMv(i)-1) = rhs_v(i);
    }
    for(int i = 0; i <= (2*m)*(2*m)-1; i++){
        rhs_v(i) = rhs_v_perm(i)*exp(COLv(i));
    }
    // End Equilibriate
    COUNT = 0;
    for(int ids = 0;ids<=2*m-1;ids++){
        for(int idr = 0;idr<=2*m-1;idr++){
            v_interp(idr,ids,0,0) = rhs_v(COUNT)*pow(hr,idr)*pow(hs,ids);
            v_interp(idr,ids,0,ns) = rhs_v(COUNT)*pow(hr,idr)*pow(hs,ids); // Use Periodicity
            COUNT++;
        }
    }

    // Right side
    rA = -0.5;
    rB = -0.5;
    sA = 0.5;
    sB = -0.5;
    // Set variable Coefficient values
    for(int idr = 0; idr <= m2p1; idr++){
        for(int ids = 0; ids <= m2p1; ids++){
            a01(idr,ids) = A01(idr,ids,nr,0)/(pow(hr,idr)*pow(hs,ids));
            a10(idr,ids) = A10(idr,ids,nr,0)/(pow(hr,idr)*pow(hs,ids));
            a11(idr,ids) = A11(idr,ids,nr,0)/(pow(hr,idr)*pow(hs,ids));
            a12(idr,ids) = A12(idr,ids,nr,0)/(pow(hr,idr)*pow(hs,ids));
            a22(idr,ids) = A22(idr,ids,nr,0)/(pow(hr,idr)*pow(hs,ids));
        }
    }
    // Boundary Conditions for u
    rhs_u.set_value(0.0);
    offset = 0;
    for(int ids = 0;ids<=m;ids++){
        for(int idr = 0;idr<=m;idr++){
            rhs_u((2*m+2)*(2*m+2)/2+offset) = tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,nr,1)/(pow(hr,idr)*pow(hs,ids));
            offset++;
            rhs_u((2*m+2)*(2*m+2)/2+offset) =  tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,nr,ns)/(pow(hr,idr)*pow(hs,ids));
            offset++;
        }
    }
    bcMat(rA,rB,sA,sB,hr,hs,a01_ptr,a10_ptr,a11_ptr,a12_ptr,a22_ptr,Mat_u_ptr,m);
    // Equilibriate
    COUNT = 0;
    for(int j = 0; j <= (2*m+2)*(2*m+2)-1;j++){
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1;i++){
            A(COUNT) = Mat_u(i,j);      
            COUNT += 1;
        }
    }
    eq_(A_ptr,&N,&NE,&LENC_number,ROW_ptr,COL_ptr,PERM_ptr);
    COUNT = 0;
    for(int j = 0; j <= (2*m+2)*(2*m+2)-1;j++){
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1;i++){
            A_mat(i,j) = Mat_u(i,j)*exp(ROW(i)+COL(j));      
        }
    }
    for(int j = 0; j <= (2*m+2)*(2*m+2)-1; j++){
        for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
            B_mat(i,j) = A_mat(i,PERM(j)-1);
        }
    }
    dgetrf_(&Mboundary,&Nboundary,B_mat_ptr,&LDAboundary,IPIV,&INFO); 
    for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
        rhs_u(i) *= exp(ROW(i));
    }
    dgetrs_(&no,&Nboundary,&ONE,B_mat_ptr,&Nboundary,IPIV,rhs_u_ptr,&Nboundary,&INFO); // Solve system
    for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
        rhs_u_perm(PERM(i)-1) = rhs_u(i);
    }
    for(int i = 0; i <= (2*m+2)*(2*m+2)-1; i++){
        rhs_u(i) = rhs_u_perm(i)*exp(COL(i));
    }
    // End Equilibriate
    COUNT = 0;
    for(int ids = 0;ids<=2*m+1;ids++){
        for(int idr = 0;idr<=2*m+1;idr++){
            u_interp(idr,ids,nr,0) = rhs_u(COUNT)*pow(hr,idr)*pow(hs,ids);
            u_interp(idr,ids,nr,ns) = rhs_u(COUNT)*pow(hr,idr)*pow(hs,ids); // Copy since periodic
            COUNT++;
        }
    }

    // Boundary Conditions for v
    rhs_v.set_value(0.0);
    offset = 0;
    for(int ids = 0; ids <= m-1; ids++){
        for(int idr = 0; idr <= m-1; idr++){
            rhs_v((2*m)*(2*m)/2+offset) = tgamma(idr+1)*tgamma(ids+1)*v(idr,ids,nr,1)/(pow(hr,idr)*pow(hs,ids));
            offset++;
            rhs_v((2*m)*(2*m)/2+offset) =  tgamma(idr+1)*tgamma(ids+1)*v(idr,ids,nr,ns)/(pow(hr,idr)*pow(hs,ids));
            offset++;
        }
    }
    bcMat(rA,rB,sA,sB,hr,hs,a01_ptr,a10_ptr,a11_ptr,a12_ptr,a22_ptr,Mat_v_ptr,m_minus);
        // Equilibriate
        COUNT = 0;
        for(int j = 0; j <= (2*m)*(2*m)-1;j++){
            for(int i = 0; i <= (2*m)*(2*m)-1;i++){
                Av(COUNT) = Mat_v(i,j);      
                COUNT += 1;
            }
        }
        eq_(Av_ptr,&Nv,&NEv,&LENC_numberv,ROWv_ptr,COLv_ptr,PERMv_ptr);
        COUNT = 0;
        for(int j = 0; j <= (2*m)*(2*m)-1;j++){
            for(int i = 0; i <= (2*m)*(2*m)-1;i++){
                A_matv(i,j) = Mat_v(i,j)*exp(ROWv(i)+COLv(j));      
            }
        }
        for(int j = 0; j <= (2*m)*(2*m)-1; j++){
            for(int i = 0; i <= (2*m)*(2*m)-1; i++){
                B_matv(i,j) = A_matv(i,PERMv(j)-1);
            }
        }
        dgetrf_(&Mboundary_v,&Nboundary_v,B_matv_ptr,&LDAboundary_v,IPIV_v,&INFO); 
        for(int i = 0; i <= (2*m)*(2*m)-1; i++){
            rhs_v(i) *= exp(ROWv(i));
        }
        dgetrs_(&no,&Nboundary_v,&ONE,B_matv_ptr,&Nboundary_v,IPIV_v,rhs_v_ptr,&Nboundary_v,&INFO); // Solve system
        for(int i = 0; i <= (2*m)*(2*m)-1; i++){
            rhs_v_perm(PERMv(i)-1) = rhs_v(i);
        }
        for(int i = 0; i <= (2*m)*(2*m)-1; i++){
            rhs_v(i) = rhs_v_perm(i)*exp(COLv(i));
        }
        // End Equilibriate
    COUNT = 0;
    for(int ids = 0;ids<=2*m-1;ids++){
        for(int idr = 0;idr<=2*m-1;idr++){
            v_interp(idr,ids,nr,0) = rhs_v(COUNT)*pow(hr,idr)*pow(hs,ids);
            v_interp(idr,ids,nr,ns) = rhs_v(COUNT)*pow(hr,idr)*pow(hs,ids); // Copy since periodic
            COUNT++;
        }
    }
}

void Hermite::BoundaryPeriodic(Darray4 &ud, Darray4 &vd){
    // LAPACK info
    int M,N,LDA,ONE;
    int Mv,Nv,LDAv;
    char no;
    double ALPHA,BETA;
    M = N = LDA = 2*m+2;
    Mv = Nv = LDAv = 2*m;
    ONE = 1;
    no = 'N';
    ALPHA = 1.0;
    BETA = 0.0;

    // Create pointers to the maps for LAPACK routines
    double *Hmap_u = Hmat_u.c_ptr();
    double *Hmap_v = Hmat_v.c_ptr();

    // Interpolation Arrays
    Darray1 uBottom, uTop;
    uBottom.define(0,2*m+1);
    uBottom.set_value(0.0);
    uTop.define(0,2*m+1);
    uTop.set_value(0.0);

    Darray1 vBottom, vTop;
    vBottom.define(0,2*m-1);
    vBottom.set_value(0.0);
    vTop.define(0,2*m-1);
    vTop.set_value(0.0);

    Darray1 uBottom_int, uTop_int;
    uBottom_int.define(0,2*m+1);
    uBottom_int.set_value(0.0);
    uTop_int.define(0,2*m+1);
    uTop_int.set_value(0.0);

    Darray1 vBottom_int, vTop_int;
    vBottom_int.define(0,2*m-1);
    vBottom_int.set_value(0.0);
    vTop_int.define(0,2*m-1);
    vTop_int.set_value(0.0);

    Darray2 uBottom_2D, uTop_2D;
    uBottom_2D.define(0,2*m+1,0,m);
    uBottom_2D.set_value(0.0);
    uTop_2D.define(0,2*m+1,0,m);
    uTop_2D.set_value(0.0);

    Darray2 vBottom_2D, vTop_2D;
    vBottom_2D.define(0,2*m-1,0,m);
    vBottom_2D.set_value(0.0);
    vTop_2D.define(0,2*m-1,0,m);
    vTop_2D.set_value(0.0);

    // Array for solution
    ud_interp.define(0,2*m+1,0,2*m+1,1,nr,1,ns);
    ud_interp.set_value(0.0);
    vd_interp.define(0,2*m-1,0,2*m-1,1,nr,1,ns);
    vd_interp.set_value(0.0);

    // Interpolate onto Boundary
    for (int i = 1; i <= nr-1; i++){
        for (int ids = 0; ids <= m; ids++){
            for (int idr = 0; idr <= m; idr++){
                uBottom(idr) = ud(idr,ids,i,ns);
                uBottom(m+1+idr) = ud(idr,ids,i+1,ns);
                uTop(idr) = ud(idr,ids,i,1);
                uTop(m+1+idr) = ud(idr,ids,i+1,1); 
            }
                // Here we interpolate onto the points (r_i+1/2,s_j) and (r_i+1/2,s_j+1)
            double *uBottom_ptr = uBottom.c_ptr();
            double *uTop_ptr = uTop.c_ptr();
            double *uBottom_int_ptr = uBottom_int.c_ptr();
            double *uTop_int_ptr = uTop_int.c_ptr();
            dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,uBottom_ptr,&ONE,&BETA,uBottom_int_ptr,&ONE);
            dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,uTop_ptr,&ONE,&BETA,uTop_int_ptr,&ONE);
            for(int i=0;i<=2*m+1;i++){
                uBottom_2D(i,ids) = uBottom_int(i);
                uTop_2D(i,ids) = uTop_int(i);
            }
        }

        // Now interpolate in y direction
        Darray1 y_interp;
        y_interp.define(0,2*m+1);
        y_interp.set_value(0.0);
        Darray1 interp;
        interp.define(0,2*m+1);
        interp.set_value(0.0);
        double *interp_ptr = interp.c_ptr();
        for(int idr = 0; idr <= 2*m+1; idr++){
            for(int ids = 0; ids <= m; ids++){
                y_interp(ids) = uBottom_2D(idr,ids);
                y_interp(m+1+ids) = uTop_2D(idr,ids);
            }
            double *y_interp_ptr = y_interp.c_ptr();
                // Here we interpolate onto the points (r_i+1/2,s_j+1/2)
            dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,y_interp_ptr,&ONE,&BETA,interp_ptr,&ONE);
            for(int ids = 0; ids <= 2*m+1; ids++){
                u_interp(idr,ids,i,0) = interp(ids);
                u_interp(idr,ids,i,ns) = interp(ids); // Copy since periodic
            }
        }
        for (int ids = 0; ids <= m-1; ids++){
            for (int idr = 0; idr <= m-1; idr++){
                vBottom(idr) = vd(idr,ids,i,ns);
                vBottom(m+idr) = vd(idr,ids,i+1,ns);
                vTop(idr) = vd(idr,ids,i,1);
                vTop(m+idr) = vd(idr,ids,i+1,1); 
            }
                // Here we interpolate onto the points (r_i+1/2,s_j) and (r_i+1/2,s_j+1)
            double *vBottom_ptr = vBottom.c_ptr();
            double *vTop_ptr = vTop.c_ptr();
            double *vBottom_int_ptr = vBottom_int.c_ptr();
            double *vTop_int_ptr = vTop_int.c_ptr();
            dgemv_(&no,&Mv,&Nv,&ALPHA,Hmap_v,&Nv,vBottom_ptr,&ONE,&BETA,vBottom_int_ptr,&ONE);
            dgemv_(&no,&Mv,&Nv,&ALPHA,Hmap_v,&Nv,vTop_ptr,&ONE,&BETA,vTop_int_ptr,&ONE);
            for(int i=0;i<=2*m-1;i++){
                vBottom_2D(i,ids) = vBottom_int(i);
                vTop_2D(i,ids) = vTop_int(i);
            }
        }

        // Now interpolate in y direction
        Darray1 y_interp_v;
        y_interp_v.define(0,2*m-1);
        y_interp_v.set_value(0.0);
        Darray1 interp_v;
        interp_v.define(0,2*m-1);
        interp_v.set_value(0.0);
        double *interp_v_ptr = interp_v.c_ptr();
        for(int idr = 0; idr <= 2*m-1; idr++){
            for(int ids = 0; ids <= m-1; ids++){
                y_interp_v(ids) = vBottom_2D(idr,ids);
                y_interp_v(m+ids) = vTop_2D(idr,ids);
            }
            double *y_interp_v_ptr = y_interp_v.c_ptr();
                // Here we interpolate onto the points (r_i+1/2,s_j+1/2)
            dgemv_(&no,&Mv,&Nv,&ALPHA,Hmap_v,&Nv,y_interp_v_ptr,&ONE,&BETA,interp_v_ptr,&ONE);
            for(int ids = 0; ids <= 2*m-1; ids++){
                v_interp(idr,ids,i,0) = interp_v(ids);
                v_interp(idr,ids,i,ns) = interp_v(ids); // Copy
            }
        }
    }
}
//--------------------------------------------------
double Bessel(double nu,double x){
    return boost::math::cyl_bessel_j(nu,x);
}

double Hermite::HBessel(double nu,double x){
    return boost::math::cyl_bessel_j(nu,x);   
}
//--------------------------------------------------
void Hermite::oversample(Darray4 &u,Darray2 &U,int nsample_r,int nsample_s){
    U.set_value(0.0);
    double rloc, sloc;
    for(int i = 1;i<=nr;i++){
        for(int j = 1;j<=ns;j++){
            for(int kr = 0;kr<=nsample_r-1;kr++){
                for(int ks = 0;ks<=nsample_s-1;ks++){
                    rloc = -0.5+double(kr)/double(nsample_r);
                    sloc = -0.5+double(ks)/double(nsample_s);
                    for(int mi = 0;mi<=2*m+1;mi++){
                        for(int mj = 0;mj<=2*m+1;mj++){
                            U(kr+(i-1)*nsample_r,ks+(j-1)*nsample_s) += u(mi,mj,i,j)*pow(rloc,mi)*pow(sloc,mj);
                        }
                    }                 
                }
            }
        }
    }

    for(int j = 1;j<=ns;j++){
        for(int ks = 0;ks<=nsample_s-1;ks++){
            sloc = -0.5+double(ks)/double(nsample_s);
            for(int mi = 0;mi<=2*m+1;mi++){
                for(int mj = 0;mj<=2*m+1;mj++){
                    U(nr*nsample_r,ks+(j-1)*nsample_s) += u(mi,mj,nr,j)*pow(0.5,mi)*pow(sloc,mj);
                }
            }
        }
    }

    for(int i = 1;i<=nr;i++){
        for(int kr =0;kr<=nsample_r-1;kr++){
            rloc = -0.5+double(kr)/double(nsample_r);
            for(int mi = 0;mi<=2*m+1;mi++){
                for(int mj = 0;mj<=2*m+1;mj++){
                    U(kr+(i-1)*nsample_r,ns*nsample_s) += u(mi,mj,i,ns)*pow(rloc,mi)*pow(0.5,mj);
                }
            }
        }
    }
    for(int mi = 0;mi<=2*m+1;mi++){
        for(int mj = 0;mj<=2*m+1;mj++){
            U(nr*nsample_r,ns*nsample_s) += u(mi,mj,nr,ns)*pow(0.5,mi)*pow(0.5,mj);
        }
    }
}


void Hermite::bcMat(double rA,double rB,double sA,double sB,double hr,double hs,
    double *a01,double* a10,double* a11,double* a12,double* a22,double *MM,int m_mat){
    if(m_mat == 1){
        matm1(&rA,&rB,&sA,&sB,a01,a10,a11,a12,a22,&hr,&hs,MM);
    }
    else if (m_mat == 2){
        matm2(&rA,&rB,&sA,&sB,a01,a10,a11,a12,a22,&hr,&hs,MM); 
    }
    else if (m_mat == 3){
        matm3(&rA,&rB,&sA,&sB,a01,a10,a11,a12,a22,&hr,&hs,MM); 
    }
}

// void Hermite::M2(double rA,double rB,double sA,double sB,double hr,double hs,
//     double *a01,double* a10,double* a11,double* a12,double* a22,double *MM){
//     matm2(&rA,&rB,&sA,&sB,a01,a10,a11,a12,a22,&hr,&hs,MM); 
// }

// void Hermite::M3(double rA,double rB,double sA,double sB,double hr,double hs,
//     double *a01,double* a10,double* a11,double* a12,double* a22,double *MM){
//     matm3(&rA,&rB,&sA,&sB,a01,a10,a11,a12,a22,&hr,&hs,MM); 
// }

