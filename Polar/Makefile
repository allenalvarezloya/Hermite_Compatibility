#
# Makefile for Hermite Boundary Conditions 
#

CXX = /usr/bin/g++ -std=c++17
FC = /usr/local/bin/gfortran
LD  = /usr/bin/g++
LDFLAGS := -lblas -llapack -lgfortran -L../BoundaryLibs/noScale/ -lHCOMPAT  
LDFLAGS += /Users/allenalvarezloya/Desktop/mc64-1.6.0/src/libmc64.a -L/usr/local/lib/gcc/x86_64-apple-darwin19/10.2.0 -L/usr/local/lib/gcc/x86_64-apple-darwin19/10.2.0/../../.. -lgfortran -lquadmath -lm 
CXXFLAGS = -O0 -DBZ_DEBUG


OBJECTS := main.o Darrays.o SetUp.o Hermite.o 
OBJECTS += polyMult.o initialData.o eqb.o allen_sub.o 
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

