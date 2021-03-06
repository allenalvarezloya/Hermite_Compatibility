SUBROUTINE get_Taylor1d(u,x,mx,mxi,hx,func,Icase)
  !
  IMPLICIT NONE
  !
  ! Uses point values evaluated by func to
  ! compute the scaled degree (mx) Taylor series
  ! centered at x
  !
  ! Here hx is the cell width and mxi is the
  ! degree of the interpolant
  !
  INTEGER, INTENT(IN) :: mx,mxi, Icase
  DOUBLE PRECISION, INTENT(IN) :: x,hx
  DOUBLE PRECISION, DIMENSION(0:mx), INTENT(OUT) :: u
  !
  INTERFACE
     DOUBLE PRECISION FUNCTION func(x,Icase)
       DOUBLE PRECISION, INTENT(IN) :: x
       INTEGER, INTENT(IN) :: Icase
     END FUNCTION func
  END INTERFACE
  !
  DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932385d0
  DOUBLE PRECISION, DIMENSION(0:mxi) :: yx,wx
  DOUBLE PRECISION :: zx
  INTEGER :: jx
  !
  ! Get data at Chebyshev nodes - note that we generate the degree mx
  ! Taylor series using degree mxi interpolants
  !
  ! Check for consistency
  !
  IF (mx > mxi) THEN
     PRINT *,'Error in get_Taylor1d - interpolation degree too low'
     STOP
  END IF

  DO jx=0,mxi
     zx=x-.5d0*hx*COS(pi*DBLE(jx)/DBLE(mxi))
     yx(jx)=func(zx,Icase)
  END DO
  !
  CALL point_to_Taylor(mxi,yx,wx)
  !
  u(:)=wx(0:mx)

END SUBROUTINE get_Taylor1d

SUBROUTINE point_to_Taylor(q,y,p)
  !
  ! Takes data on q+1 Chebyshev nodes and
  ! uses the Bjorck-Peyreyra algorithm to express
  ! the interpolant as a degree q Taylor series
  ! (see Dahlquist-Bjorck volume 1 page 375)
  !
  INTEGER, INTENT(IN) :: q
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(IN) :: y
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(OUT) :: p
  !
  DOUBLE PRECISION, DIMENSION(0:q) :: z
  DOUBLE PRECISION :: pi=3.1415926535897932385d0
  INTEGER :: j,k
  !
  DO j=0,q
     z(j)=-.5d0*COS(pi*DBLE(j)/DBLE(q))
  END DO
  !
  ! First divided differences
  !
  p=y
  DO k=1,q
     DO j=q,k,-1
        p(j)=(p(j)-p(j-1))/(z(j)-z(j-k))
     END DO
  END DO
  !
  ! Now Horner
  !
  DO k=q-1,0,-1
     DO j=k,q-1
        p(j)=p(j)-z(k)*p(j+1)
     END DO
  END DO
  !
END SUBROUTINE point_to_Taylor

SUBROUTINE indat(floc,r,hr,m,Icase)
  implicit none
  INTEGER, INTENT(IN) :: m, Icase
  DOUBLE PRECISION, INTENT(IN) :: hr,r
  DOUBLE PRECISION, DIMENSION(0:m) :: floc
  
  INTERFACE 
     DOUBLE PRECISION FUNCTION ss(r,Icase)
       DOUBLE PRECISION, INTENT(IN) :: r 
       INTEGER, INTENT(IN) :: Icase
     END FUNCTION ss
  END INTERFACE
  
  CALL get_Taylor1d(floc,r,m,20,hr,ss,Icase)

END SUBROUTINE indat

DOUBLE PRECISION FUNCTION ss(r,Icase)
  !
  DOUBLE PRECISION, INTENT(IN) :: r
  DOUBLE PRECISION :: x,y 
  INTEGER, INTENT(IN) :: Icase
  DOUBLE PRECISION :: pi = 3.1415926535897932385d0
  ! Cartesian
  ! x = r 
  ! y = s
  ! Polar Grid
  x = 2*r-1
  IF (Icase == 0) THEN
     ss = exp(-500.d0*(x+0.5d0)**2)
     ! ss = sin(pi*x)
  ELSE IF (Icase == 1) THEN
     ss = 0.0 !1000.d0*(x+0.5d0)*exp(-500.d0*(x+0.5d0)**2)
  ELSE if (Icase == 2) THEN
     ss = x 
  END IF
  !
END FUNCTION ss
