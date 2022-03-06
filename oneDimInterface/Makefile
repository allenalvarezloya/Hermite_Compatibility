#
# Makefile for Hermite Boundary Conditions 
#

CXX = g++
FC  = gfortran -ffree-line-length-none
LD  = g++
LDFLAGS = -lgfortran -lblas -llapack 
FCFLAGS = -O0 -fno-range-check
CXXINCLUDE = 
CXXFLAGS = -O0 
AR = ar 

FILE = obj1.txt
OBJECTS_M1 :=$(file < $(FILE))

FILE = obj2.txt
OBJECTS_M2 :=$(file < $(FILE))

FILE = obj3.txt
OBJECTS_M3 :=$(file < $(FILE))

FILE = intObj1.txt
OBJECTS_intM1 :=$(file < $(FILE))

FILE = intObj2.txt
OBJECTS_intM2 :=$(file < $(FILE))

FILE = intObj3.txt
OBJECTS_intM3 :=$(file < $(FILE))

# FILE = obj4.txt
# OBJECTS_M4 :=$(file < $(FILE))

# FILE = obj5.txt
# OBJECTS_M5 :=$(file < $(FILE))


EXEC = main.x	

.PHONY: clean

all: libM1 libM2 libM3 libint1 libint2 libint3 #libM4 libM5

libM1: $(OBJECTS_M1) 
	$(AR) rcs libHCOMPAT.a $(OBJECTS_M1) 

libM2: $(OBJECTS_M2) 
	$(AR) rcs libHCOMPAT.a $(OBJECTS_M2) 


libM3: $(OBJECTS_M3) 
	$(AR) rcs libHCOMPAT.a $(OBJECTS_M3) 

libint1: $(OBJECTS_intM1) 
	$(AR) rcs libHCOMPAT.a $(OBJECTS_intM1) 

libint2: $(OBJECTS_intM2) 
	$(AR) rcs libHCOMPAT.a $(OBJECTS_intM2) 

libint3: $(OBJECTS_intM3) 
	$(AR) rcs libHCOMPAT.a $(OBJECTS_intM3) 

# libM4: $(OBJECTS_M4) 
# 	$(AR) rcs libHCOMPAT.a $(OBJECTS_M4) 

# libM5: $(OBJECTS_M5) 
# 	$(AR) rcs libHCOMPAT.a $(OBJECTS_M5) 


exe: $(EXEC)


$(EXEC): $(OBJECTS)
	$(LD) $(LDFLAGS) $(OBJECTS) -o $(EXEC)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -c $<

%.o : %.f90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

clean:
	rm -f $(OBJECTS_M1) $(EXEC) 

sandwich: 
	make -j 15

