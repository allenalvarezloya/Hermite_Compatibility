SUBROUTINE tpolymul(a,b,c,q)
  !
  !  truncated multiplication of degree q polynomials
  !  with coefficients a,b. On return c will contain the
  !  coefficients of the product truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  INTEGER :: j,k
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(IN) :: a,b
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(OUT) :: c
  !
  c=0.d0
  DO j=0,q
     DO k=0,j
        c(j)=c(j)+a(k)*b(j-k)
     END DO
  END DO

END SUBROUTINE tpolymul

SUBROUTINE tpolymul2(a,b,c,q)
  !
  !  truncated multiplication of tensor product
  !  degree q polynomials in 2 variables
  !  with coefficients a,b. On return c will contain the
  !  coefficients of the product truncated at degree q
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: q
  DOUBLE PRECISION, DIMENSION(0:q,0:q), INTENT(IN) :: a,b
  DOUBLE PRECISION, DIMENSION(0:q,0:q), INTENT(OUT) :: c
  INTEGER :: j1,j2,k1,k2
  !
  c=0.d0
  DO j2=0,q
     DO j1=0,q
        DO k2=0,j2
           DO k1=0,j1
              c(j1,j2)=c(j1,j2)+a(k1,k2)*b(j1-k1,j2-k2)
           END DO
        END DO
     END DO
  END DO
END SUBROUTINE tpolymul2

SUBROUTINE tpolypow2(a,c,r,q)
 !
 !  truncated power function of a tensor product degree
 !  (q,q) polynomial
 !  with coefficients a. On return c will contain the
 !  coefficients of a**r truncated at degree q
 !
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: q
 DOUBLE PRECISION, DIMENSION(0:q,0:q), INTENT(IN) :: a
 DOUBLE PRECISION, DIMENSION(0:q,0:q), INTENT(OUT) :: c
 DOUBLE PRECISION, INTENT(IN) :: r
 INTEGER :: j1,j2,k1,k2
 !
 c=0.d0
 c(0,0)=a(0,0)**r
 DO k1=1,q
  DO k2=1,k1
   c(0,k1)=c(0,k1)+a(0,k2)*c(0,k1-k2)*(DBLE(k2)*(1.d0+r)-DBLE(k1))
 END DO
 c(0,k1)=c(0,k1)/(DBLE(k1)*a(0,0))
END DO
DO j1=1,q
  DO j2=1,j1
   c(j1,0)=c(j1,0)+a(j2,0)*c(j1-j2,0)*(DBLE(j2)*(1.d0+r)-DBLE(j1))
 END DO
 c(j1,0)=c(j1,0)/(DBLE(j1)*a(0,0))
 DO k1=1,q
   DO j2=1,j1
    DO k2=0,k1
     c(j1,k1)=c(j1,k1)+   &
     a(j2,k2)*c(j1-j2,k1-k2)*(DBLE(j2)*(1.d0+r)-DBLE(j1))
   END DO
 END DO
 DO k2=1,k1
  c(j1,k1)=c(j1,k1)-DBLE(j1)*a(0,k2)*c(j1,k1-k2)
END DO
c(j1,k1)=c(j1,k1)/(DBLE(j1)*a(0,0))
END DO
END DO

END SUBROUTINE tpolypow2