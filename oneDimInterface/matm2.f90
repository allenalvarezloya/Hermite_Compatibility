subroutine matm2(rA,hr,M2) bind(c) 
implicit none 
 double precision :: rA, hr 
 double precision :: M2(0:5,0:5) 
 M2(0,0) = 1 
 M2(0,1) = 0 
 M2(0,2) = 0 
 M2(0,3) = 0 
 M2(0,4) = 0 
 M2(0,5) = 0 
 M2(1,0) = 0 
 M2(1,1) = 0 
 M2(1,2) = 2 
 M2(1,3) = 0 
 M2(1,4) = 0 
 M2(1,5) = 0 
 M2(2,0) = 0 
 M2(2,1) = 0 
 M2(2,2) = 0 
 M2(2,3) = 0 
 M2(2,4) = 24 
 M2(2,5) = 0 
 M2(3,0) = 1 
 M2(3,1) = rA*hr 
 M2(3,2) = hr**2*rA**2 
 M2(3,3) = hr**3*rA**3 
 M2(3,4) = hr**4*rA**4 
 M2(3,5) = hr**5*rA**5 
 M2(4,0) = 0 
 M2(4,1) = 1 
 M2(4,2) = 2*rA*hr 
 M2(4,3) = 3*hr**2*rA**2 
 M2(4,4) = 4*hr**3*rA**3 
 M2(4,5) = 5*hr**4*rA**4 
 M2(5,0) = 0 
 M2(5,1) = 0 
 M2(5,2) = 2 
 M2(5,3) = 6*rA*hr 
 M2(5,4) = 12*hr**2*rA**2 
 M2(5,5) = 20*hr**3*rA**3 
end subroutine