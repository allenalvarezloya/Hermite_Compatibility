SUBROUTINE get_Taylor2d(u,x,y,mx,my,mxi,myi,hx,hy,func,Icase,t)
  !
  IMPLICIT NONE
  !
  ! Uses point values evaluated by func to
  ! compute the scaled degree (mx,my) Taylor series
  ! centered at x,y
  !
  ! Here hx X hy is the cell width and mxi,myi are the
  ! degrees of the interpolant
  !
  INTEGER, INTENT(IN) :: mx,my,mxi,myi, Icase
  DOUBLE PRECISION, INTENT(IN) :: x,hx,y,hy,t
  DOUBLE PRECISION, DIMENSION(0:mx,0:my), INTENT(OUT) :: u
  !
  INTERFACE
     DOUBLE PRECISION FUNCTION func(x,y,Icase,t)
       DOUBLE PRECISION, INTENT(IN) :: x,y,t
       INTEGER, INTENT(IN) :: Icase
     END FUNCTION func
  END INTERFACE
  !
  DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932385d0
  DOUBLE PRECISION, DIMENSION(0:mxi) :: yx,wx
  DOUBLE PRECISION, DIMENSION(0:myi) :: yy,wy
  DOUBLE PRECISION, DIMENSION(0:mx,0:myi) :: ux
  DOUBLE PRECISION :: zx,zy
  INTEGER :: jx,jy
  !
  ! Get data at Chebyshev nodes - note that we generate the degree mx,my
  ! Taylor series using degree mxi,myi interpolants
  !
  ! Check for consistency
  !
  IF ((mx > mxi).OR.(my > myi)) THEN
     PRINT *,'Error in get_Taylor2d - interpolation degree too low'
     STOP
  END IF

  DO jy=0,myi
     zy=y-.5d0*hy*COS(pi*DBLE(jy)/DBLE(myi))
     DO jx=0,mxi
        zx=x-.5d0*hx*COS(pi*DBLE(jx)/DBLE(mxi))
        yx(jx)=func(zx,zy,Icase,t)
     END DO
     !
     CALL point_to_Taylor(mxi,yx,wx)
     !
     ux(:,jy)=wx(0:mx)
  END DO
  !
  DO jx=0,mx
     yy=ux(jx,:)
     CALL point_to_Taylor(myi,yy,wy)
     u(jx,:)=wy(0:my)
  END DO

END SUBROUTINE get_Taylor2d

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

SUBROUTINE indat(floc,r,s,hr,hs,m,Icase,t)

  implicit none
  INTEGER, INTENT(IN) :: m, Icase
  DOUBLE PRECISION, INTENT(IN) :: hs,hr,r,s,t
  DOUBLE PRECISION, DIMENSION(0:m,0:m) :: floc
  
  INTERFACE 
     DOUBLE PRECISION FUNCTION ss(r,s,Icase,t)
       DOUBLE PRECISION, INTENT(IN) :: r,s,t 
       INTEGER, INTENT(IN) :: Icase
     END FUNCTION ss
  END INTERFACE
  
  CALL get_Taylor2d(floc,r,s,m,m,15,15,hr,hs,ss,Icase,t)

END SUBROUTINE indat

DOUBLE PRECISION FUNCTION ss(r,s,Icase,t)
  !
  DOUBLE PRECISION, INTENT(IN) :: r,s,t
  DOUBLE PRECISION :: x,y 
  INTEGER, INTENT(IN) :: Icase
  DOUBLE PRECISION :: pi = 3.1415926535897932385d0
  DOUBLE PRECISION :: k31 = 6.3801618959239383506237d0;
  DOUBLE PRECISION :: k33 = 13.01520072169843441983d0;
  DOUBLE PRECISION :: gamma = 6.3801618959239383506237d0/13.01520072169843441983d0;
  DOUBLE PRECISION :: alpha = 1d0-6.3801618959239383506237d0/13.01520072169843441983d0;
  ! Cartesian
  ! x = r 
  ! y = s
  ! Polar Grid
  x = (alpha*r + gamma)*cos(2d0*pi*s)
  y = (alpha*r + gamma)*sin(2d0*pi*s)
  IF (Icase == 0) THEN
   ss = BESSEL_JN(3, k33*(alpha*r+gamma))*cos(6d0*pi*s)*cos(k33*t) 
   ! ss = sin(2d0*pi*x)*sin(2d0*pi*y)*cos(2d0**1.5*pi*t)
  ELSE IF (Icase == 1) THEN
   ss = x 
  ELSE IF (Icase == 2) THEN
   ss = y 
  END IF 
  !
END FUNCTION ss