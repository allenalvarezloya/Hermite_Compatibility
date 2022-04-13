subroutine intm0(hu,hv,C,D,M0) bind(c) 
implicit none 
 double precision :: hu,hv,C,D 
 double precision :: M0(0:3,0:3) 
 M0(0,0) = 1 
 M0(0,1) = 0 
 M0(0,2) = -1 
 M0(0,3) = 0 
 M0(1,0) = 0 
 M0(1,1) = C**2/hu 
 M0(1,2) = 0 
 M0(1,3) = -D**2/hv 
 M0(2,0) = 1 
 M0(2,1) = -.5 
 M0(2,2) = 0 
 M0(2,3) = 0 
 M0(3,0) = 0 
 M0(3,1) = 0 
 M0(3,2) = 1 
 M0(3,3) = .5 
end subroutine