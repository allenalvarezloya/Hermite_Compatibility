#include "Hermite.h"
using namespace std;
//------------------------------------------
//------------------------------------------
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
//------------------------------------------
//------------------------------------------
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
//------------------------------------------
//------------------------------------------
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
//------------------------------------------ 
void Hermite::interpPD(Darray2 &u){
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

    // Arrays for interpolation
    Darray1 U;
    U.define(0,2*m+1);
    U.set_value(0.0);
    double *U_ptr = U.c_ptr();
    ud_interp.define(0,2*m+1,1,nr);
    ud_interp.set_value(0.0);
    for(int i = 0; i <= nr-1; i++){
        for(int idr = 0; idr <= m; idr++){
            U(idr) = u(idr,i);
            U(m+1+idr) = u(idr,i+1);
        }
        dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,U_ptr,&ONE,&BETA,U_ptr,&ONE);
        for(int idr= 0; idr <= 2*m+1; idr++){
            ud_interp(idr,i+1) = U(idr);
        }  
    }
}
//------------------------------------------
//------------------------------------------ 
void Hermite::interpDP(Darray2 &ud){
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

    // Arrays for interpolation
    Darray1 U;
    U.define(0,2*m+1);
    U.set_value(0.0);
    double *U_ptr = U.c_ptr();

    u_interp.define(0,2*m+1,0,nr);
    u_interp.set_value(0.0);
    for(int i = 1; i <= nr-1; i++){
        for(int idr = 0; idr <= m; idr++){
            U(idr) = ud(idr,i);
            U(m+1+idr) = ud(idr,i+1);
        }
        dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,U_ptr,&ONE,&BETA,U_ptr,&ONE);
        for(int idr = 0; idr <= 2*m+1; idr++){
            u_interp(idr,i) = U(idr);
        }  
    }
}
//------------------------------------------
//------------------------------------------ 
void Hermite::boundaryReflection(Darray2 &ud){
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

    // Arrays for interpolation
    Darray1 U;
    U.define(0,2*m+1);
    U.set_value(0.0);
    double *U_ptr = U.c_ptr();

    for(int idr = 0; idr <= m; idr++){
        U(idr) = pow(-1.0,idr+1)*ud(idr,1);
        U(m+1+idr) = ud(idr,1);
    }
    dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,U_ptr,&ONE,&BETA,U_ptr,&ONE);
    for(int idr = 0; idr <= 2*m+1; idr++){
        u_interp(idr,0) = U(idr);
    }  
    for(int idr = 0; idr <= m; idr++){
        U(idr) = ud(idr,nr);
        U(m+1+idr) = pow(-1.0,idr+1)*ud(idr,nr);
    }
    dgemv_(&no,&M,&N,&ALPHA,Hmap_u,&N,U_ptr,&ONE,&BETA,U_ptr,&ONE);
    for(int idr = 0; idr <= 2*m+1; idr++){
        u_interp(idr,nr) = U(idr);
    }  
}
//------------------------------------------
//------------------------------------------ 
void Hermite::laplacian(int m,Darray1 &u,double h){
    int m2p1 = 2*m+1;
    Darray1 lapU;
    lapU.define(0,m2p1);
    lapU.set_value(0.0);
    for(int idr = 0; idr <= m2p1-2; idr++){
        lapU(idr) = double(idr+2)*double(idr+1)*u(idr+2)/(h*h); 
    }
    u.set_value(0.0);
    for(int idr = 0; idr <= m2p1-2; idr++){
        u(idr) = lapU(idr); 
    }
}
//------------------------------------------
//------------------------------------------ 
void Hermite::recursion(Darray2 &u_int, Darray2 &um_int,Darray2 &u,int i,double h,double C,double dt){
    int m2p1 = 2*m+1;
    Darray1 Sum;
    Sum.define(0,m2p1);
    Darray1 lapU;
    lapU.define(0,m2p1);
    Sum.set_value(0.0);
    lapU.set_value(0.0);
    for(int idr = 0; idr <= m2p1; idr++){
        Sum(idr) = u_int(idr,i);
        lapU(idr) = u_int(idr,i);
    }
    laplacian(m,lapU,h);
    double scale;
    for(int l = 1; l <= 2*m; l++){
        scale = pow(0.5*C*dt,2.0*l)/tgamma(2.0*l+1);
        for(int idr = 0; idr <= m2p1; idr++){
            Sum(idr) += scale*lapU(idr);
        }
        laplacian(m,lapU,h);
    }
    for(int idr = 0; idr <= m; idr++){
        u(idr,i) = 2*Sum(idr) - um_int(idr,i);
    }
}
//------------------------------------------
//------------------------------------------ 
void Hermite::recursion2(Darray2 &u_int,Darray2 &u,int i,double h,double C,double dt){
    double S = 0.0;
    double scale = 0.0;
    for(int k = 0; k <= m; k++){
        S = 0.0;
        for(int l = 0; l <= m - floor(k/2); l++){
            scale = pow(0.5*C*(dt/h),2*l)*tgamma(2*l+k+1)/tgamma(2*l+1)/tgamma(k+1);
            S += scale*u_int(2*l+k,i);
        }
        u(k,i) = 2*S - u(k,i);
    }
}
//------------------------------------------
//------------------------------------------ 
void Hermite::boundaryConditions(int m, double hr,Darray2 &u,double t){
    int m2p1 = 2*m+1;
    // LAPACK info for u
    int Mboundary,Nboundary,LDAboundary;
    Mboundary = Nboundary = LDAboundary = (2*m+2);
    int IPIV[(2*m+2)];
    int INFO;
    int ONE = 1;
    char no = 'N';
    char fName[100];

    Darray2 Mat_u;
    Mat_u.define(0,m2p1,0,m2p1);
    Mat_u.set_value(0.0);
    double* Mat_u_ptr = Mat_u.c_ptr();

    Darray1 rhs_u;
    rhs_u.define(0,m2p1);
    double* rhs_u_ptr = rhs_u.c_ptr();

    // Left side
    double rA = 0.5;
    rhs_u.set_value(0.0);
    boundaryValues(m,hr,t,rhs_u);
    for(int idr = 0; idr <= m; idr++){
        rhs_u(m+1+idr) = tgamma(idr+1)*u(idr,1);
    }
    bcMat(rA,hr,Mat_u_ptr,m);
    sprintf(fName, "M.ext");
    outPutMat(Mat_u,0,2*m+1,0,2*m+1,fName);
    dgetrf_(&Mboundary,&Nboundary,Mat_u_ptr,&LDAboundary,IPIV,&INFO);
    dgetrs_(&no,&Nboundary,&ONE,Mat_u_ptr,&Nboundary,IPIV,rhs_u_ptr,&Nboundary,&INFO);
    for(int idr = 0; idr <= m2p1; idr++){
        u_interp(idr,0) = rhs_u(idr);
    }
    // Right side
    rA = -0.5;
    rhs_u.set_value(0.0);
    boundaryValues(m,hr,t,rhs_u);
    for(int idr = 0; idr <= m; idr++){
        rhs_u(m+1+idr) = tgamma(idr+1)*u(idr,nr);
    }
    bcMat(rA,hr,Mat_u_ptr,m);
    dgetrf_(&Mboundary,&Nboundary,Mat_u_ptr,&LDAboundary,IPIV,&INFO);
    dgetrs_(&no,&Nboundary,&ONE,Mat_u_ptr,&Nboundary,IPIV,rhs_u_ptr,&Nboundary,&INFO);
    for(int idr = 0; idr <= m2p1; idr++){
        u_interp(idr,nr) = rhs_u(idr);
    }
}
//------------------------------------------
//------------------------------------------ 
void Hermite::boundaryValues(int m,double h,double t,Darray1 &RHS){
    double pi = 3.1415926535897932385;
    for(int idr = 0; idr <= m; idr ++){
        RHS(idr) = pow(-1.0,idr)*pow(2*pi*h,2*idr)*sin(-2*pi*t);
    }
}
//------------------------------------------
//------------------------------------------ 
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
//------------------------------------------
//------------------------------------------ 
void Hermite::outPutVec(Darray1 &M,int i_start,int i_end,char* fName){
    FILE *extFile = fopen(fName, "w");
    for (int i=i_start; i<=i_end; i++){
            fprintf(extFile, "%18.10e \n ", M(i));
        }
    fclose(extFile);
}
//------------------------------------------
//------------------------------------------ 
void Hermite::bcMat(double rA,double hr,double *MM,int m_mat){
    if(m_mat == 1){
        matm1(&rA,&hr,MM);
    }
    else if (m_mat == 2){
        matm2(&rA,&hr,MM);
    }
    else if (m_mat == 3){
        matm3(&rA,&hr,MM);
    }
}
//------------------------------------------
//------------------------------------------ 
// void Hermite::interface(double hu,double hv,double C,double D,double *MM, int m_mat){
//     if(m_mat == 1){
//         intm1(&hu,&hv,&C,&D,MM);
//     }
//     else if(m_mat == 2){
//         intm2(&hu,&hv,&C,&D,MM);
//     }
//     else if(m_mat == 3){
//         intm3(&hu,&hv,&C,&D,MM);
//     }
// }

void Hermite::oversample(Darray2 &u,int nr,int m,char* fName){
    int nrefine = 10;
    double h_loc = 1.0/double(nrefine);
    double x_loc;
    double sum = 0.0;
    int count = 0;
    Darray1 u_oversample;
    u_oversample.define(0,nr*nrefine);
    u_oversample.set_value(0.0);
    for (int i = 1; i <= nr; i++){
        for(int j = 0; j <= nrefine-1; j++){
            x_loc = -0.5 + h_loc*double(j);
            for(int idr = 0; idr <= 2*m+1; idr ++){
                sum += u(idr,i)*pow(x_loc,idr);
            }
            u_oversample(count) = sum;
            count += 1;
            sum = 0.0;
        }
    }
    outPutVec(u_oversample,0,nr*nrefine,fName);
}

void Hermite::arrayCopy(Darray2 &u,Darray2 &v,int is,int ie,int js,int je){
    for(int i = is; i <= ie; i++){
        for(int j = js; j <= je; j++){
            u(i,j) = v(i,j);
        }
    }
}

