SUBROUTINE get_Taylor1d(u,x,t,mx,mxi,hx,func)
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
  INTEGER, INTENT(IN) :: mx,mxi
  DOUBLE PRECISION, INTENT(IN) :: x,hx,t
  DOUBLE PRECISION, DIMENSION(0:mx), INTENT(OUT) :: u
  !
  INTERFACE
     DOUBLE PRECISION FUNCTION func(x,t)
       DOUBLE PRECISION, INTENT(IN) :: x, t
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
     yx(jx)=func(zx,t)
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

SUBROUTINE indat(floc,r,t,hr,m)
  implicit none
  INTEGER, INTENT(IN) :: m
  DOUBLE PRECISION, INTENT(IN) :: hr,r,t
  DOUBLE PRECISION, DIMENSION(0:m) :: floc
  
  INTERFACE 
     DOUBLE PRECISION FUNCTION ss(r,t)
       DOUBLE PRECISION, INTENT(IN) :: r, t
     END FUNCTION ss
  END INTERFACE
  
  CALL get_Taylor1d(floc,r,t,m,15,hr,ss)

END SUBROUTINE indat

DOUBLE PRECISION FUNCTION ss(r,t)
  !
  DOUBLE PRECISION, INTENT(IN) :: r, t
  DOUBLE PRECISION :: x,y 
  DOUBLE PRECISION :: pi = 3.1415926535897932385d0
  ! Cartesian
  ! x = r 
  ! y = s
  ! Polar Grid
  x = r
  ss = sin(2*pi*x)*cos(2*pi*t)
  !
END FUNCTION ss