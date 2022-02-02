#pragma once
/* The actual coefficients of the differencing on the mesh,
	and the influence of the spatial variation of epsilon:  */
#include "MHD.h"

void init_solve();
void tridag(double* a, double* b, double* c, double* r, double* utri, double* gam, int n);
int solve(double u_in[ZMAX][RMAX], double s[ZMAX][RMAX], int itermax, double tol_test);


void potential_boundary();
void electric_field();
void potential_solve();

void mag_phi();