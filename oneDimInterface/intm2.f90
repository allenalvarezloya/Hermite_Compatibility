subroutine intm2(hu,hv,C,D,M2) bind(c) 
implicit none 
 double precision :: hu,hv,C,D 
 double precision :: M2(0:11,0:11) 
 M2(0,0) = 1 
 M2(0,1) = 0 
 M2(0,2) = 0 
 M2(0,3) = 0 
 M2(0,4) = 0 
 M2(0,5) = 0 
 M2(0,6) = -1 
 M2(0,7) = 0 
 M2(0,8) = 0 
 M2(0,9) = 0 
 M2(0,10) = 0 
 M2(0,11) = 0 
 M2(1,0) = 0 
 M2(1,1) = C**2/hu 
 M2(1,2) = 0 
 M2(1,3) = 0 
 M2(1,4) = 0 
 M2(1,5) = 0 
 M2(1,6) = 0 
 M2(1,7) = -D**2/hv 
 M2(1,8) = 0 
 M2(1,9) = 0 
 M2(1,10) = 0 
 M2(1,11) = 0 
 M2(2,0) = 0 
 M2(2,1) = 0 
 M2(2,2) = 2/hu**2 
 M2(2,3) = 0 
 M2(2,4) = 0 
 M2(2,5) = 0 
 M2(2,6) = 0 
 M2(2,7) = 0 
 M2(2,8) = -2/hv**2 
 M2(2,9) = 0 
 M2(2,10) = 0 
 M2(2,11) = 0 
 M2(3,0) = 0 
 M2(3,1) = 0 
 M2(3,2) = 0 
 M2(3,3) = 6*C**2/hu**3 
 M2(3,4) = 0 
 M2(3,5) = 0 
 M2(3,6) = 0 
 M2(3,7) = 0 
 M2(3,8) = 0 
 M2(3,9) = -6*D**2/hv**3 
 M2(3,10) = 0 
 M2(3,11) = 0 
 M2(4,0) = 0 
 M2(4,1) = 0 
 M2(4,2) = 0 
 M2(4,3) = 0 
 M2(4,4) = 24/hu**4 
 M2(4,5) = 0 
 M2(4,6) = 0 
 M2(4,7) = 0 
 M2(4,8) = 0 
 M2(4,9) = 0 
 M2(4,10) = -24/hv**4 
 M2(4,11) = 0 
 M2(5,0) = 0 
 M2(5,1) = 0 
 M2(5,2) = 0 
 M2(5,3) = 0 
 M2(5,4) = 0 
 M2(5,5) = 120*C**2/hu**5 
 M2(5,6) = 0 
 M2(5,7) = 0 
 M2(5,8) = 0 
 M2(5,9) = 0 
 M2(5,10) = 0 
 M2(5,11) = -120*D**2/hv**5 
 M2(6,0) = 1 
 M2(6,1) = -.5 
 M2(6,2) = .25 
 M2(6,3) = -.125 
 M2(6,4) = .625e-1 
 M2(6,5) = -.3125e-1 
 M2(6,6) = 0 
 M2(6,7) = 0 
 M2(6,8) = 0 
 M2(6,9) = 0 
 M2(6,10) = 0 
 M2(6,11) = 0 
 M2(7,0) = 0 
 M2(7,1) = 1/hu 
 M2(7,2) = -1.0/hu 
 M2(7,3) = .75/hu 
 M2(7,4) = -.500/hu 
 M2(7,5) = .3125/hu 
 M2(7,6) = 0 
 M2(7,7) = 0 
 M2(7,8) = 0 
 M2(7,9) = 0 
 M2(7,10) = 0 
 M2(7,11) = 0 
 M2(8,0) = 0 
 M2(8,1) = 0 
 M2(8,2) = 2/hu**2 
 M2(8,3) = -3.0/hu**2 
 M2(8,4) = 3.00/hu**2 
 M2(8,5) = -2.500/hu**2 
 M2(8,6) = 0 
 M2(8,7) = 0 
 M2(8,8) = 0 
 M2(8,9) = 0 
 M2(8,10) = 0 
 M2(8,11) = 0 
 M2(9,0) = 0 
 M2(9,1) = 0 
 M2(9,2) = 0 
 M2(9,3) = 0 
 M2(9,4) = 0 
 M2(9,5) = 0 
 M2(9,6) = 1 
 M2(9,7) = .5 
 M2(9,8) = .25 
 M2(9,9) = .125 
 M2(9,10) = .625e-1 
 M2(9,11) = .3125e-1 
 M2(10,0) = 0 
 M2(10,1) = 0 
 M2(10,2) = 0 
 M2(10,3) = 0 
 M2(10,4) = 0 
 M2(10,5) = 0 
 M2(10,6) = 0 
 M2(10,7) = 1/hv 
 M2(10,8) = 1.0/hv 
 M2(10,9) = .75/hv 
 M2(10,10) = .500/hv 
 M2(10,11) = .3125/hv 
 M2(11,0) = 0 
 M2(11,1) = 0 
 M2(11,2) = 0 
 M2(11,3) = 0 
 M2(11,4) = 0 
 M2(11,5) = 0 
 M2(11,6) = 0 
 M2(11,7) = 0 
 M2(11,8) = 2/hv**2 
 M2(11,9) = 3.0/hv**2 
 M2(11,10) = 3.00/hv**2 
 M2(11,11) = 2.500/hv**2 
end subroutine