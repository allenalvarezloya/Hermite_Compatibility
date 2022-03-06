#
# Makefile for Hermite Boundary Conditions 
#

CXX = g++ -std=c++17 -I/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk/usr/include
FC = gfortran
LD  = g++
LDFLAGS := -lgfortran -lblas -llapack -L../ -lHCOMPAT
CXXFLAGS = -O0 -DBZ_DEBUG


OBJECTS := main.o Darrays.o SetUp.o initialData.o Hermite.o
EXEC = main.x	

.PHONY: clean

exe: $(EXEC)


$(EXEC): $(OBJECTS)
	$(LD) $(LDFLAGS) $(OBJECTS) -o $(EXEC)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDE) -c $<

%.o : %.f90
	$(FC) $(FCFLAGS) $(FCINCLUDE) -c $<

clean:
	rm -f $(OBJECTS) $(EXEC) 

sandwich: 
	make clean
	make 
	./main.x

