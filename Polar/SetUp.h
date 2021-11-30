#ifndef __SETUP_H_INCLUDED__

#define __SETUP_H_INCLUDED__
#include <cmath>
#include <iostream>
#include "Darrays.h"
#include <boost/math/special_functions/bessel.hpp>
using namespace std;
extern "C" void tpolymul2_(double *,double *,double *, int *);
extern "C" void tpolypow2_(double *,double *,double *, int *);
class SetUp {
public:
 	int m;                         // Number of derivatives used
 	int nr, ns;                    // Number of cells in r and s
    int metricOrder;               // Order we need for the metric      
 	double hr, hs;                 // Stepsize in r and s
 	double x_s, x_e, y_s, y_e;     // Start and endpoints 
 	void compute_grids();          // Compute grids in r and s
 	void transformation();         // Computes the transformation
    void transformation_noScale(); // Computes the transformation (no scaling)
 	void compute_initial_data(Darray4 &u,Darray4 &v); //initial data
    void compute_initial_data_dual_ext(Darray4 &u,Darray4 &v); //initial data
 	void Metric_Info();            // Compute metric info
    void Metric_Info_noScale();            // Compute metric info
    void metricCoeffs();           // Compute Metric coeffs
    void metricCoeffs_noScale();           // Compute Metric coeffs
 	Darray1 r,s,rd,sd;             // Allocate r and s grids
    Darray1 rd_ext, sd_ext;        // Allocate r and s extended grids
 	Darray1 x,y;                   // Allocate x and y grids
    Darray4 X,Y,Xd,Yd;             // Allocate physical coordinates 
    Darray4 X_ns, Y_ns;             // Allocate X no scaling Y no scaling
    Darray4 rx,ry,sx,sy;           // Allocate metric info (primal)
    Darray4 rx_ns, ry_ns, sx_ns, sy_ns; // No scaling
    Darray4 rxd,ryd,sxd,syd;       // Allocate metric info (dual)
    Darray4 rxx,ryy,sxx,syy;       // Allocate Metric info (second ders primal)
    Darray4 a01, a10, a11,a12,a22; // Variable Coefficients
    Darray4 a01_ns, a10_ns, 
    a11_ns,a12_ns,a22_ns;          // Variable Coefficients
    Darray4 a01d, a10d, 
    a11d,a12d,a22d;                // Variable Coefficients on dual
};

#endif