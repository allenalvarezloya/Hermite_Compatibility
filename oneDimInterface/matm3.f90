subroutine matm3(rA,hr,M3) bind(c) 
implicit none 
 double precision :: rA, hr 
 double precision :: M3(0:7,0:7) 
 M3(0,0) = 1 
 M3(0,1) = 0 
 M3(0,2) = 0 
 M3(0,3) = 0 
 M3(0,4) = 0 
 M3(0,5) = 0 
 M3(0,6) = 0 
 M3(0,7) = 0 
 M3(1,0) = 0 
 M3(1,1) = 0 
 M3(1,2) = 2 
 M3(1,3) = 0 
 M3(1,4) = 0 
 M3(1,5) = 0 
 M3(1,6) = 0 
 M3(1,7) = 0 
 M3(2,0) = 0 
 M3(2,1) = 0 
 M3(2,2) = 0 
 M3(2,3) = 0 
 M3(2,4) = 24 
 M3(2,5) = 0 
 M3(2,6) = 0 
 M3(2,7) = 0 
 M3(3,0) = 0 
 M3(3,1) = 0 
 M3(3,2) = 0 
 M3(3,3) = 0 
 M3(3,4) = 0 
 M3(3,5) = 0 
 M3(3,6) = 720 
 M3(3,7) = 0 
 M3(4,0) = 1 
 M3(4,1) = rA 
 M3(4,2) = rA**2 
 M3(4,3) = rA**3 
 M3(4,4) = rA**4 
 M3(4,5) = rA**5 
 M3(4,6) = rA**6 
 M3(4,7) = rA**7 
 M3(5,0) = 0 
 M3(5,1) = 1 
 M3(5,2) = 2*rA 
 M3(5,3) = 3*rA**2 
 M3(5,4) = 4*rA**3 
 M3(5,5) = 5*rA**4 
 M3(5,6) = 6*rA**5 
 M3(5,7) = 7*rA**6 
 M3(6,0) = 0 
 M3(6,1) = 0 
 M3(6,2) = 2 
 M3(6,3) = 6*rA 
 M3(6,4) = 12*rA**2 
 M3(6,5) = 20*rA**3 
 M3(6,6) = 30*rA**4 
 M3(6,7) = 42*rA**5 
 M3(7,0) = 0 
 M3(7,1) = 0 
 M3(7,2) = 0 
 M3(7,3) = 6 
 M3(7,4) = 24*rA 
 M3(7,5) = 60*rA**2 
 M3(7,6) = 120*rA**3 
 M3(7,7) = 210*rA**4 
end subroutine