#ifndef __SETUP_H_INCLUDED__

#define __SETUP_H_INCLUDED__
#include <cmath>
#include <iostream>
#include "Darrays.h"

using namespace std;
extern "C" void indat_(double *,double *,double *,double *,int *);

class SetUp{
public:
	int m;                       // Number of derivatives
	int nr;                      // Number of cells
	double hr;                   // Cell width
	Darray1 r,rd;                // Allocate grid arrays
	void compute_grids();        // Compute primal and dual grids
	void compute_initial_data(Darray2 &u,Darray1 &r,int ns,int ne,double t); // Compute iniital data
};

#endif