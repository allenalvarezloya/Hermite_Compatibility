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
//--------------------------------------------------
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
void Hermite::interpolate_pd(Darray4 &u){
    // This subroutine interpolates from primal to dual

    // LAPACK info
    int M,N,LDA,ONE;
    char no;
    double ALPHA,BETA;
    M = N = LDA = 2*m+2;
    ONE = 1;
    no = 'N';
    ALPHA = 1.0;
    BETA = 0.0;

    // Create pointers to the maps for LAPACK routines
    double *Hmap_u = Hmat_u.c_ptr();

    // Interpolation Arrays
    Darray1 uBottom, uTop;
    uBottom.define(0,2*m+1);
    uBottom.set_value(0.0);
    uTop.define(0,2*m+1);
    uTop.set_value(0.0);

    Darray1 uBottom_int, uTop_int;
    uBottom_int.define(0,2*m+1);
    uBottom_int.set_value(0.0);
    uTop_int.define(0,2*m+1);
    uTop_int.set_value(0.0);

    Darray2 uBottom_2D, uTop_2D;
    uBottom_2D.define(0,2*m+1,0,m);
    uBottom_2D.set_value(0.0);
    uTop_2D.define(0,2*m+1,0,m);
    uTop_2D.set_value(0.0);

    // Array for solution
    ud_interp.define(0,2*m+1,0,2*m+1,1,nr,1,ns);
    ud_interp.set_value(0.0);

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
        }
    }
}
//--------------------------------------------------
void Hermite::interpolate_dp(Darray4 &ud){
    // This subroutine interpolates from dual to priaml

    // LAPACK info
    int M,N,LDA,ONE;
    char no;
    double ALPHA,BETA;
    M = N = LDA = 2*m+2;
    ONE = 1;
    no = 'N';
    ALPHA = 1.0;
    BETA = 0.0;

    // Create pointers to the maps for LAPACK routines
    double *Hmap_u = Hmat_u.c_ptr();

    // Interpolation Arrays
    Darray1 uBottom, uTop;
    uBottom.define(0,2*m+1);
    uBottom.set_value(0.0);
    uTop.define(0,2*m+1);
    uTop.set_value(0.0);

    Darray1 uBottom_int, uTop_int;
    uBottom_int.define(0,2*m+1);
    uBottom_int.set_value(0.0);
    uTop_int.define(0,2*m+1);
    uTop_int.set_value(0.0);

    Darray2 uBottom_2D, uTop_2D;
    uBottom_2D.define(0,2*m+1,0,m);
    uBottom_2D.set_value(0.0);
    uTop_2D.define(0,2*m+1,0,m);
    uTop_2D.set_value(0.0);

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
void Hermite::recursion_curvilinear(Darray4 &u_int,int i,int j,
    Darray4 &A01, Darray4 &A10,Darray4 & A11,Darray4 &A12,Darray4 & A22,Darray4 &u,
    double hr,double hs,double dt){
    int m2p1 = 2*m+1;
    // Define arrays needed for computation
    Darray2 lapSum;
    lapSum.define(0,m2p1,0,m2p1);
    lapSum.set_value(0.0);
    double *lapSum_ptr = lapSum.c_ptr();
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

    for(int ids = 0; ids <= m2p1; ids++){
        for(int idr = 0; idr <= m2p1; idr++){
            lapSum(idr,ids) = u_int(idr,ids,i,j);
        }
    }
    Darray2 u_loc;
    u_loc.define(0,m2p1,0,m2p1);
    for(int ids = 0; ids <= m2p1; ids++){
        for(int idr = 0; idr <= m2p1; idr++){
            u_loc(idr,ids) = u_int(idr,ids,i,j);
        }
    }
    curviLaplacian(hr,hs,m,a01_recursion,a10_recursion,a11_recursion,a12_recursion,a22_recursion,u_loc);
    for(int k = 1; k <= 2*m; k++){
        for(int ids = 0; ids <= m2p1; ids++){
            for(int idr = 0; idr <= m2p1; idr++){
                lapSum(idr,ids) += pow(0.5*dt,2*k)*u_loc(idr,ids)/tgamma(2*k+1);
            }
        }
        curviLaplacian(hr,hs,m,a01_recursion,a10_recursion,a11_recursion,a12_recursion,a22_recursion,u_loc);
    }
    for(int ids = 0; ids <= m; ids++){
        for(int idr = 0; idr <= m; idr++){
            u(idr,ids,i,j) = 2*lapSum(idr,ids) - u(idr,ids,i,j);
        }
    }
}
//--------------------------------------------------
void Hermite::curviLaplacian(double hx,double hy,int m,Darray2 &a01,Darray2 &a10,Darray2 &a11,Darray2 &a12,Darray2 &a22,Darray2 &u){
    int m2p1 = 2*m+1;
    Darray2 uxx,uyy,ux,uy,uxy;
    uxx.define(0,m2p1,0,m2p1);
    uyy.define(0,m2p1,0,m2p1);
    ux.define(0,m2p1,0,m2p1);
    uy.define(0,m2p1,0,m2p1);
    uxy.define(0,m2p1,0,m2p1);
    uxx.set_value(0.0);
    uyy.set_value(0.0);
    ux.set_value(0.0);
    uy.set_value(0.0);
    uxy.set_value(0.0);

    double* uxx_ptr = uxx.c_ptr();
    double* uyy_ptr = uyy.c_ptr();
    double* ux_ptr = ux.c_ptr();
    double* uy_ptr = uy.c_ptr();
    double* uxy_ptr = uxy.c_ptr();

    double* a11_ptr = a11.c_ptr();
    double* a22_ptr = a22.c_ptr();
    double* a10_ptr = a10.c_ptr();
    double* a01_ptr = a01.c_ptr();
    double* a12_ptr = a12.c_ptr();

    Darray2 a11_uxx, a22_uyy, a12_uxy, a10_ux, a01_uy;
    a11_uxx.define(0,m2p1,0,m2p1);
    a22_uyy.define(0,m2p1,0,m2p1);
    a10_ux.define(0,m2p1,0,m2p1);
    a01_uy.define(0,m2p1,0,m2p1);
    a12_uxy.define(0,m2p1,0,m2p1);
    a11_uxx.set_value(0.0);
    a22_uyy.set_value(0.0);
    a10_ux.set_value(0.0);
    a01_uy.set_value(0.0);
    a12_uxy.set_value(0.0);


    double* a11_uxx_ptr = a11_uxx.c_ptr();
    double* a22_uyy_ptr = a22_uyy.c_ptr();
    double* a10_ux_ptr = a10_ux.c_ptr();
    double* a01_uy_ptr = a01_uy.c_ptr();
    double* a12_uxy_ptr = a12_uxy.c_ptr();
    // uxx
    for(int j = 0; j <= m2p1; j++){
        for(int i = 0; i <= m2p1 - 2; i++){
            uxx(i,j) = (i+1)*(i+2)*u(i+2,j)/(hx*hx);
        }
    }
    // uyy
    for(int i = 0; i <= m2p1; i++){
        for(int j = 0; j <= m2p1 - 2; j++){
            uyy(i,j) = (j+1)*(j+2)*u(i,j+2)/(hy*hy);
        }
    }
    // ux 
    for(int j = 0; j <= m2p1; j++){
        for(int i = 0; i <= m2p1 - 1; i++){
            ux(i,j) = (i+1)*u(i+1,j)/hx;
        }
    }
    // uy
    for(int i = 0; i <= m2p1; i++){
        for(int j = 0; j <= m2p1 - 1; j++){
            uy(i,j) = (j+1)*u(i,j+1)/hy;
        }
    }
    // uxy
    for(int i = 0; i <= m2p1 - 1; i++){
        for(int j = 0; j <= m2p1 - 1; j++){
            uxy(i,j) = (i+1)*(j+1)*u(i+1,j+1)/(hx*hy);
        }
    }
    tpolymul2_(a11_ptr,uxx_ptr,a11_uxx_ptr,&m2p1);
    tpolymul2_(a22_ptr,uyy_ptr,a22_uyy_ptr,&m2p1);
    tpolymul2_(a10_ptr,ux_ptr,a10_ux_ptr,&m2p1);
    tpolymul2_(a01_ptr,uy_ptr,a01_uy_ptr,&m2p1);
    tpolymul2_(a12_ptr,uxy_ptr,a12_uxy_ptr,&m2p1);
    // Compute Laplacian
    for(int j = 0; j <= m2p1; j++){
        for(int i = 0; i <= m2p1; i++){
            u(i,j) = a11_uxx(i,j) + a22_uyy(i,j) + a10_ux(i,j) + a01_uy(i,j) + a12_uxy(i,j);
        }
    }
}
//--------------------------------------------------
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
//--------------------------------------------------
void Hermite::outPutVec(Darray1 &M,int i_start,int i_end,char* fName){
    FILE *extFile = fopen(fName, "w");
    for (int i=i_start; i<=i_end; i++){
            fprintf(extFile, "%18.10e \n ", M(i));
        }
    fclose(extFile);
}
//--------------------------------------------------
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
    }
}
//--------------------------------------------------
void Hermite::boundaryConditions(int m,double hr,double hs,Darray4 &u, 
    Darray4 &A01,Darray4 & A10,Darray4 & A11,Darray4 &A12,Darray4 & A22){
    int m2p1 = 2*m+1;
    // LAPACK info for u
    int Mboundary,Nboundary,LDAboundary;
    Mboundary = Nboundary = LDAboundary = (2*m+2)*(2*m+2);
    int IPIV[(2*m+2)*(2*m+2)];
    int INFO;
    int ONE = 1;
    char no = 'N';
    char fName[100];

    // Define arrays
    Darray2 a10, a01, a11, a12, a22;
    a10.define(0,1,0,1);
    a01.define(0,1,0,1);
    a11.define(0,1,0,1);
    a12.define(0,1,0,1);
    a22.define(0,1,0,1);
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
    // Declare parameters
    double rA,rB,sA,sB;
    Darray1 rhs_u;
    rhs_u.define(0,(2*m+2)*(2*m+2)-1);
    double* rhs_u_ptr = rhs_u.c_ptr();
    int offset = 0;
    int COUNT = 0;

    Darray1 rhs_u_perm;
    rhs_u_perm.define(0,(2*m+2)*(2*m+2)-1);
    rhs_u_perm.set_value(0.0);

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

    for(int j = 1; j <= ns - 1; j++){
        // Left side
        rA = 0.5;
        rB = 0.5;
        sA = 0.5;
        sB = -0.5;
        // Set variable Coefficient values
        for(int idr = 0; idr <= 1; idr++){
            for(int ids = 0; ids <= 1; ids++){
                a01(idr,ids) = A01(idr,ids,0,j);
                a10(idr,ids) = A10(idr,ids,0,j);
                a11(idr,ids) = A11(idr,ids,0,j);
                a12(idr,ids) = A12(idr,ids,0,j);
                a22(idr,ids) = A22(idr,ids,0,j);
            }
        }
        // Boundary Conditions for u
        rhs_u.set_value(0.0);
        offset = 0;
        for(int ids = 0;ids<=m;ids++){
            for(int idr = 0;idr<=m;idr++){
                rhs_u((2*m+2)*(2*m+2)/2+offset) = tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,1,j+1);
                offset++;
                rhs_u((2*m+2)*(2*m+2)/2+offset) = tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,1,j);
                offset++;
            }
        }
        bcMat(rA,rB,sA,sB,hr,hs,a01_ptr,a10_ptr,a11_ptr,a12_ptr,a22_ptr,Mat_u_ptr,m);
        sprintf(fName, "M%0.5d.ext",j);
        outPutMat(Mat_u,0,(2*m+2)*(2*m+2)-1,0,(2*m+2)*(2*m+2)-1,fName);
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
        // dgetrf_(&Mboundary,&Nboundary,Mat_u_ptr,&LDAboundary,IPIV,&INFO); 
        // dgetrs_(&no,&Nboundary,&ONE,Mat_u_ptr,&Nboundary,IPIV,rhs_u_ptr,&Nboundary,&INFO); // Solve system
        COUNT = 0;
        for(int ids = 0;ids<=2*m+1;ids++){
            for(int idr = 0;idr<=2*m+1;idr++){
                u_interp(idr,ids,0,j) = rhs_u(COUNT);
                COUNT++;
            }
        }





        // Right side
        rA = -0.5;
        rB = -0.5;
        sA = 0.5;
        sB = -0.5;
        // Set variable Coefficient values
        for(int idr = 0; idr <= 1; idr++){
            for(int ids = 0; ids <= 1; ids++){
                a01(idr,ids) = A01(idr,ids,nr,j);
                a10(idr,ids) = A10(idr,ids,nr,j);
                a11(idr,ids) = A11(idr,ids,nr,j);
                a12(idr,ids) = A12(idr,ids,nr,j);
                a22(idr,ids) = A22(idr,ids,nr,j);
            }
        }
        // Boundary Conditions for u
        rhs_u.set_value(0.0);
        offset = 0;
        for(int ids = 0;ids<=m;ids++){
            for(int idr = 0;idr<=m;idr++){
                rhs_u((2*m+2)*(2*m+2)/2+offset) = tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,nr,j+1);
                offset++;
                rhs_u((2*m+2)*(2*m+2)/2+offset) =  tgamma(idr+1)*tgamma(ids+1)*u(idr,ids,nr,j);
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
        // dgetrf_(&Mboundary,&Nboundary,Mat_u_ptr,&LDAboundary,IPIV,&INFO); 
        // dgetrs_(&no,&Nboundary,&ONE,Mat_u_ptr,&Nboundary,IPIV,rhs_u_ptr,&Nboundary,&INFO); // Solve system
        COUNT = 0;
        for(int ids = 0;ids<=2*m+1;ids++){
            for(int idr = 0;idr<=2*m+1;idr++){
                u_interp(idr,ids,nr,j) = rhs_u(COUNT);
                COUNT++;
            }
        }
    }
}
//--------------------------------------------------
void Hermite::BoundaryPeriodic(Darray4 &ud){
    // LAPACK info
    int M,N,LDA,ONE;
    char no;
    double ALPHA,BETA;
    M = N = LDA = 2*m+2;
    ONE = 1;
    no = 'N';
    ALPHA = 1.0;
    BETA = 0.0;

    // Create pointers to the maps for LAPACK routines
    double *Hmap_u = Hmat_u.c_ptr();

    // Interpolation Arrays
    Darray1 uBottom, uTop;
    uBottom.define(0,2*m+1);
    uBottom.set_value(0.0);
    uTop.define(0,2*m+1);
    uTop.set_value(0.0);

    Darray1 uBottom_int, uTop_int;
    uBottom_int.define(0,2*m+1);
    uBottom_int.set_value(0.0);
    uTop_int.define(0,2*m+1);
    uTop_int.set_value(0.0);

    Darray2 uBottom_2D, uTop_2D;
    uBottom_2D.define(0,2*m+1,0,m);
    uBottom_2D.set_value(0.0);
    uTop_2D.define(0,2*m+1,0,m);
    uTop_2D.set_value(0.0);

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
    }
    // Start Corners 
    // Interpolate onto Boundary
    for (int ids = 0; ids <= m; ids++){
        for (int idr = 0; idr <= m; idr++){
            uBottom(idr) = ud(idr,ids,nr,ns);
            uBottom(m+1+idr) = ud(idr,ids,1,ns);
            uTop(idr) = ud(idr,ids,nr,1);
            uTop(m+1+idr) = ud(idr,ids,1,1); 
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
            u_interp(idr,ids,0,0) = interp(ids);
            u_interp(idr,ids,0,ns) = interp(ids);  // Copy since periodic
            u_interp(idr,ids,nr,0) = interp(ids); // Copy since periodic
            u_interp(idr,ids,nr,ns) = interp(ids);  // Copy since periodic
        }
    }
    // End Corners
}
//--------------------------------------------------
double Bessel(double nu,double x){
    return boost::math::cyl_bessel_j(nu,x);
}
//--------------------------------------------------
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
//--------------------------------------------------
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
//--------------------------------------------------