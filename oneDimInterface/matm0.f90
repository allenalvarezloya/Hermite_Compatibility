subroutine matm0(rA,hr,M0) bind(c) 
implicit none 
 double precision :: rA, hr 
 double precision :: M0(0:1,0:1) 
 M0(0,0) = 1 
 M0(0,1) = 0 
 M0(1,0) = 1 
 M0(1,1) = rA 
end subroutine