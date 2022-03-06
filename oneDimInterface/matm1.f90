subroutine matm1(rA,hr,M1) bind(c) 
implicit none 
 double precision :: rA, hr 
 double precision :: M1(0:3,0:3) 
 M1(0,0) = 1 
 M1(0,1) = 0 
 M1(0,2) = 0 
 M1(0,3) = 0 
 M1(1,0) = 0 
 M1(1,1) = 0 
 M1(1,2) = 2 
 M1(1,3) = 0 
 M1(2,0) = 1 
 M1(2,1) = rA*hr 
 M1(2,2) = hr**2*rA**2 
 M1(2,3) = hr**3*rA**3 
 M1(3,0) = 0 
 M1(3,1) = 1 
 M1(3,2) = 2*rA*hr 
 M1(3,3) = 3*hr**2*rA**2 
end subroutine