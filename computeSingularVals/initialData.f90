SUBROUTINE get_Taylor2d(u,x,y,mx,my,mxi,myi,hx,hy,func,Icase)
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
  DOUBLE PRECISION, INTENT(IN) :: x,hx,y,hy
  DOUBLE PRECISION, DIMENSION(0:mx,0:my), INTENT(OUT) :: u
  !
  INTERFACE
     DOUBLE PRECISION FUNCTION func(x,y,Icase)
       DOUBLE PRECISION, INTENT(IN) :: x,y
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
        yx(jx)=func(zx,zy,Icase)
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

SUBROUTINE indat(floc,r,s,hr,hs,m,Icase)
  ! This routine outputs:
  !  rho(x,y,t) = \sin(x) \sin(y) \sin(\omega t), with \omega = \sqrt(2)
  !  u(x,y,t) = 1/\omega \cos(x) \sin(y) \cos(\omega t),
  !  v(x,y,t) = 1/\omega \sin(x) \cos(y) \cos(\omega t).

  implicit none
  INTEGER, INTENT(IN) :: m, Icase
  DOUBLE PRECISION, INTENT(IN) :: hs,hr,r,s
  DOUBLE PRECISION, DIMENSION(0:m,0:m) :: floc
  
  INTERFACE 
     DOUBLE PRECISION FUNCTION ss(r,s,Icase)
       DOUBLE PRECISION, INTENT(IN) :: r,s 
       INTEGER, INTENT(IN) :: Icase
     END FUNCTION ss
  END INTERFACE
  
  CALL get_Taylor2d(floc,r,s,m,m,2*m+1,2*m+1,hr,hs,ss,Icase)

END SUBROUTINE indat

DOUBLE PRECISION FUNCTION ss(r,s,Icase)
  !
  DOUBLE PRECISION, INTENT(IN) :: r,s
  DOUBLE PRECISION :: x,y 
  INTEGER, INTENT(IN) :: Icase
  DOUBLE PRECISION :: pi = 3.1415926535897932385d0
  DOUBLE PRECISION :: k31 = 6.3801618959239383506237d0;
  DOUBLE PRECISION :: k33 = 13.01520072169843441983d0;
  DOUBLE PRECISION :: gamma = 6.3801618959239383506237d0/13.01520072169843441983d0;
  DOUBLE PRECISION :: alpha = 1-6.3801618959239383506237d0/13.01520072169843441983d0;
  ! Cartesian
  x = r 
  y = s
  ! Polar Grid
  ! x = (alpha*r + gamma)*cos(2*pi*s)
  ! y = (alpha*r + gamma)*sin(2*pi*s)
  IF (Icase == 0) THEN
   ss = BESSEL_JN(3, k33*(alpha*r+gamma))*cos(6*pi*s) 
  ELSE IF (Icase == 1) THEN
   ss = x 
  ELSE IF (Icase == 2) THEN
   ss = y 
  END IF 
  !
END FUNCTION ss