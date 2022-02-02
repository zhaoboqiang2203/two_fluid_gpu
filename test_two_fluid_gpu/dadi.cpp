
/*
====================================================================

DADIXY.CPP

  This is a dynamic ADI solution of the equation
  Del( epsilon * Del(phi) ) = rho
  in z-r coordinates.
  Neumann, Dirichlet, and symmetry boundary conditions are supported.
  Symmetry boundary condition is only applied when r0=0, and the
  neuman flag is specified there with zero slope on the electric field.

  dxdxu + dydyu =  s

  The function is based on a Peaceman Rachford Douglas
  advance of the parabolic equation:

  dtu = dxdxu + dydyu -s

  But with the time step dynamically chosen so as to speed convergence
  to the dtu=0 steady state solution which is what we want.  It is the
  user's responsiblity to put the initial guess of the solution stored
  in u when the function is called.  u=0 or u=u(previous time step)
  are possible choices.

  The function sends back the finishing iteration number, the finishing
  normalized maximum error and location, and the failure code.

  The failure codes are:
	 ifail=0, success actually
	 ifail=1, had to discard too many times
	 ifail=2, used up the iterations

  Send in a negative tol_test if you want to freeze the time step to
  the initial choice.  Send in a 0 adel_t to let program pick an
  initial del_t.

Revision/Programmer/Date
0.?	(Peterm, ??-??-94)	Conversion from C DADI.
0.9	(JohnV Peterm 08-09-95) Modify epsilon in set_coeff for dielectric
		blocks.

====================================================================
*/
#include "dadi.h"
#include <math.h>


double a_x1geom[ZMAX][RMAX], b_x1geom[ZMAX][RMAX], c_x1geom[ZMAX][RMAX];
double a_x2geom[ZMAX][RMAX], b_x2geom[ZMAX][RMAX], c_x2geom[ZMAX][RMAX];

/*  The arrays used internal to DADI which contain the coefficients
	  for each tridiagonal matrix solution */
double a_x1[ZMAX], b_x1[ZMAX], c_x1[ZMAX];
double a_x2[ZMAX], b_x2[ZMAX], c_x2[ZMAX];
/*  Various copies of the 'answer' we're working toward */
double u[ZMAX][RMAX];
double uwork[ZMAX][RMAX], ustor[ZMAX][RMAX], ustar[ZMAX][RMAX];
double r_x1[ZMAX], v_x1[ZMAX], gam_x1[ZMAX];
double r_x2[ZMAX], v_x2[ZMAX], gam_x2[ZMAX];
/*  Our fictitious time step */
double del_t0;
/*  epsilon at the grid locations */
//double** epsi;
/*  The size of the system */

#ifndef MAX
#define MAX(x, y)       (((x) > (y)) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y)       (((x) < (y)) ? (x) : (y))
#endif

#ifndef DBL_MIN
#define DBL_MIN         1E-200
#endif

double Er[ZMAX][RMAX];
double Ez[ZMAX][RMAX];



/**********************************************************************
  Single Peaceman Rachford Douglas pass with Direchlet 0 c boundary
  conditions for the equation:

  dtu = dxdxu + dydyu -s, where s is constant in time.

  The Crank-Nicolson finite difference approximation to the
  above equation leads to the fractional step or
  ADI equations:

  u*(i,j)-(del_t/2dxdx)[u*(i+1,j)-2u*(i,j)+u*(i-1,j)]
  = un(i,j)+(del_t/2dydy)[un(i,j+1)-2un(i,j)+un(i,j-1)] - (del_t/2)s(i,j)

  un+1(i,j)-(del_t/2dydy)[un+1(i,j+1)-2un+1(i,j)+un+1(i,j-1)]
  = u*(i,j)+(del_t/2dxdx)[u*(i+1,j)-2u*(i,j)+u*(i-1,j)] - (del_t/2)s(i,j)

  **********************************************************************/



  /******************************************************/

void adi(double uadi[ZMAX][RMAX], double s[ZMAX][RMAX], double del_t)
{
	register int i, j;
	double dth;


	dth = .5 * del_t;

	/***************************************/
	/* Do z pass.  Set up variables for    */
	/* tridiagonal inversion along z.      */

	for (j = 0; j < RMAX; j++)
	{
		for (i = 0; i < ZMAX; i++)
		{
			a_x1[i] = -dth * a_x1geom[i][j];
			b_x1[i] = 1 - dth * b_x1geom[i][j];
			c_x1[i] = -dth * c_x1geom[i][j];
		}

		/*  handle the boundary conditions, neumann and dirichlet */
		for (i = 0; i < ZMAX; i++)
			if (b_x2geom[i][j] != 0)  //non-dirichlet
				r_x1[i] = uadi[i][j] + dth * (-s[i][j]
					+ ((j > 0) ? a_x2geom[i][j] * uadi[i][j - 1] : 0)
					+ b_x2geom[i][j] * uadi[i][j]
					+ ((j < RMAX) ? c_x2geom[i][j] * uadi[i][j + 1] : 0));
			else
				r_x1[i] = uadi[i][j];  // metal sets the potential here.

	  /* Solve tridiagonal system. */
		tridag(a_x1, b_x1, c_x1, r_x1, v_x1, gam_x1, ZMAX);

		/* Copy solution into ustar. */
		for (i = 0; i < ZMAX; i++) ustar[i][j] = v_x1[i];
	}

	/***************************************/
	/* Do y pass.  Set up variables for    */
	/* tridiagonal inversion along y     */

	for (i = 0; i < ZMAX; i++)
	{
		for (j = 0; j < RMAX; j++)
		{
			a_x2[j] = -dth * a_x2geom[i][j];
			b_x2[j] = 1 - dth * b_x2geom[i][j];
			c_x2[j] = -dth * c_x2geom[i][j];
		}

		/*  handle the boundary conditions, dirichlet or neumann */

		 /*  The following code handles some special cases like corners*/
		for (j = 0; j < RMAX; j++)
			if (b_x2geom[i][j] != 0)  //non-dirichlet
				r_x2[j] = ustar[i][j] + dth * (-s[i][j]
					+ ((i > 0) ? a_x1geom[i][j] * ustar[i - 1][j] : 0)
					+ b_x1geom[i][j] * ustar[i][j]
					+ ((i < ZMAX) ? c_x1geom[i][j] * ustar[i + 1][j] : 0));
			else
				r_x2[j] = ustar[i][j];  // metal sets the potential here.

	  /* Solve tridiagonal system. */
		tridag(a_x2, b_x2, c_x2, r_x2, v_x2, gam_x2, RMAX);

		/* Copy solution into ustar. */
		for (j = 0; j < RMAX; j++) uadi[i][j] = v_x2[j];

	}

	/*****************************************************/
	/* Dirchlet boundary conditions for i=0 and i=RMAX. */

}
void init_solve()
{
	register int i, j, k;

	del_t0 = dt;
	for (i = 0; i < ZMAX; i++)
	{
		for (j = 0; j < RMAX; j++)
		{
			if (btype[i][j] == 1)
			{
				a_x1geom[i][j] = 1 / (dz * dz);
				c_x1geom[i][j] = 1 / (dz * dz);
				b_x1geom[i][j] = -(a_x1geom[i][j] + c_x1geom[i][j]);

				a_x2geom[i][j] = 1 / (dr * dr) - 1 / (2 * j * dr * dr);
				c_x2geom[i][j] = 1 / (dr * dr) + 1 / (2 * j * dr * dr);
				b_x2geom[i][j] = -(a_x2geom[i][j] + c_x2geom[i][j]);
			}
			else if (btype[i][j] == 110)
			{
				if (i == 0)
				{
					a_x1geom[i][j] = 0;
				}
				else
				{
					a_x1geom[i][j] = 1 / (dz * dz);
				}
				
				if (i == ZMAX - 1)
				{
					c_x1geom[i][j] = 0;
				}
				else
				{
					c_x1geom[i][j] = 1 / (dz * dz);
				}

				b_x1geom[i][j] = -(a_x1geom[i][j] + c_x1geom[i][j]);

				if (j == 0)
				{
					a_x2geom[i][j] = 0;
					c_x2geom[i][j] = 1 / (dr * dr) + 1 / (dr * dr);
					b_x2geom[i][j] = -(a_x2geom[i][j] + c_x2geom[i][j]);
				}
				else if (j == RMAX - 1)
				{
					a_x2geom[i][j] = 1 / (dr * dr) - 1 / (2 * j * dr * dr);
					c_x2geom[i][j] = 0;
					b_x2geom[i][j] = -(a_x2geom[i][j] + c_x2geom[i][j]);
				}
				else
				{
					a_x2geom[i][j] = 1 / (dr * dr) + 1 / (2 * j * dr * dr);
					c_x2geom[i][j] = 1 / (dr * dr) - 1 / (2 * j * dr * dr);
					b_x2geom[i][j] = -(a_x2geom[i][j] + c_x2geom[i][j]);
				}
				
			}
			else if(ptype[i][j] == CATHODE_BOUNDARY || ptype[i][j] == INLET)
			{
				if (i == 0)
				{
					a_x1geom[i][j] = 0;
				}
				else
				{
					a_x1geom[i][j] = 1 / (dz * dz);
				}

				if (i == ZMAX - 1)
				{
					c_x1geom[i][j] = 0;
				}
				else
				{
					c_x1geom[i][j] = 1 / (dz * dz);
				}

				b_x1geom[i][j] = -(a_x1geom[i][j] + c_x1geom[i][j]);

				if (j == 0)
				{
					a_x2geom[i][j] = 0;
					c_x2geom[i][j] = 1 / (dr * dr) + 1 / (dr * dr);
					b_x2geom[i][j] = -(a_x2geom[i][j] + c_x2geom[i][j]);
				}
				else if (j == RMAX - 1)
				{
					a_x2geom[i][j] = 1 / (dr * dr) - 1 / (2 * j * dr * dr);
					c_x2geom[i][j] = 0;
					b_x2geom[i][j] = -(a_x2geom[i][j] + c_x2geom[i][j]);
				}
				else
				{
					a_x2geom[i][j] = 1 / (dr * dr) - 1 / (2 * j * dr * dr);
					c_x2geom[i][j] = 1 / (dr * dr) + 1 / (2 * j * dr * dr);
					b_x2geom[i][j] = -(a_x2geom[i][j] + c_x2geom[i][j]);
				}

			}
			else if (ptype[i][j] == ANODE_BOUNDARY)
			{
				a_x1geom[i][j] = 0.0;
				b_x1geom[i][j] = 0.0;
				c_x1geom[i][j] = 0.0;
				a_x2geom[i][j] = 0.0;
				b_x2geom[i][j] = 0.0;
				c_x2geom[i][j] = 0.0;
			}
			else if (ptype[i][j] == CYLINDRICAL_AXIS)
			{
				a_x1geom[i][j] = 1 / (dz * dz);
				c_x1geom[i][j] = 1 / (dz * dz);
				b_x1geom[i][j] = -(a_x1geom[i][j] + c_x1geom[i][j]);

				a_x2geom[i][j] = 0;
				c_x2geom[i][j] = 1 / (dr * dr);
				b_x2geom[i][j] = -(a_x2geom[i][j] + c_x2geom[i][j]);

			}
			else
			{
				if (i == 0)
				{
					a_x1geom[i][j] = 0;
				}
				else
				{
					a_x1geom[i][j] = 1 / (dz * dz);
				}

				if (i == ZMAX - 1)
				{
					c_x1geom[i][j] = 0;
				}
				else
				{
					c_x1geom[i][j] = 1 / (dz * dz);
				}

				b_x1geom[i][j] = -(a_x1geom[i][j] + c_x1geom[i][j]);


				if (j == 0)
				{
					a_x2geom[i][j] = 0;
					c_x2geom[i][j] = 1 / (dr * dr) + 1 / (dr * dr);
					b_x2geom[i][j] = -(a_x2geom[i][j] + c_x2geom[i][j]);
				}
				else if (j == RMAX - 1)
				{
					a_x2geom[i][j] = 1 / (dr * dr) - 1 / (2 * j * dr * dr);
					c_x2geom[i][j] = 0;
					b_x2geom[i][j] = -(a_x2geom[i][j] + c_x2geom[i][j]);
				}
				else
				{
					a_x2geom[i][j] = 1 / (dr * dr) - 1 / (2 * j * dr * dr);
					c_x2geom[i][j] = 1 / (dr * dr) + 1 / (2 * j * dr * dr);
					b_x2geom[i][j] = -(a_x2geom[i][j] + c_x2geom[i][j]);
				}
			}

		}
	}

	
#ifdef DADI_DEBUG

	matrix_to_csv((double**)a_x1geom, ZMAX, RMAX, RMAX, (char*)(".\\output\\a_x1geom.csv"));
	matrix_to_csv((double**)b_x1geom, ZMAX, RMAX, RMAX, (char*)(".\\output\\b_x1geom.csv"));
	matrix_to_csv((double**)c_x1geom, ZMAX, RMAX, RMAX, (char*)(".\\output\\c_x1geom.csv"));
	matrix_to_csv((double**)a_x2geom, ZMAX, RMAX, RMAX, (char*)(".\\output\\a_x2geom.csv"));
	matrix_to_csv((double**)b_x2geom, ZMAX, RMAX, RMAX, (char*)(".\\output\\b_x2geom.csv"));
	matrix_to_csv((double**)c_x2geom, ZMAX, RMAX, RMAX, (char*)(".\\output\\c_x2geom.csv"));
#endif
}


/**********************************************************************/

/*dadi(u_in, s, itermax, tol_test, u_x0, u_xlx, u_y0, u_yly)
int itermax;
double tol_test;
double **u_in, **s, *u_x0, *u_xlx, *u_y0, *u_yly;
*/
int solve(double u_in[ZMAX][RMAX], double s[ZMAX][RMAX], int itermax, double tol_test)
{
	register int i, j;
	int iter, ndiscard;
	static double del_t = 0.0;
	double del_td = 0, tptop = 0, tpbot = 0, ratio = 0;
	double rnorm = 0, rsum = 0, res = 0, errchk = 0, dxdxutrm = 0, dydyutrm = 0;

	rnorm = rsum = 0.0;
	for (i = 0; i < ZMAX; i++)
		for (j = 0; j < RMAX; j++) {

			/* Residual normalization.  */
			// dirichlet conditions don't add to rnorm
			rnorm += ((b_x2geom[i][j] == 0) ? 0 : s[i][j] * s[i][j]);

			/*  copy u_in to u for working purposes.  */
			u[i][j] = (double)u_in[i][j];

			//calculate an initial estimate of the residual
		/* Use the residual as the absolute error and if it is bigger
		   than the stored maximum absolute error, record new maximum
		   absolute error and location.  */

			if (i > 0 && j > 0 && i < ZMAX && j < RMAX) {
				dxdxutrm = a_x1geom[i][j] * u_in[i - 1][j] + b_x1geom[i][j] * u_in[i][j] + c_x1geom[i][j] * u_in[i + 1][j];
				dydyutrm = a_x2geom[i][j] * u_in[i][j - 1] + b_x2geom[i][j] * u_in[i][j] + c_x2geom[i][j] * u_in[i][j + 1];
			}

			/* only include points which are not in structures. */
			errchk = ((b_x2geom[i][j] == 0) ? 0 : dxdxutrm + dydyutrm - s[i][j]);

			/* Residual sums. */
			rsum += errchk * errchk;
		}

	// If rnorm is zero, we must deal with it...
	if (rnorm == 0.0) {

		// check dirichlet conditions
		for (i = 0; i < ZMAX; i++) for (j = 0; j < RMAX; j++)
		{
			rnorm += sqr(((i > 0 && b_x2geom[i - 1][j] != 0) ? c_x1geom[i - 1][j] * u[i][j] : 0));
			// check right
			rnorm += sqr(((i < (ZMAX - 1) && b_x2geom[i + 1][j] != 0) ? a_x1geom[i + 1][j] * u[i][j] : 0));
			// check up
			rnorm += sqr(((j > 0 && b_x2geom[i][j - 1] != 0) ? c_x2geom[i][j - 1] * u[i][j] : 0));
			// check right
			rnorm += sqr(((j < (RMAX - 1) && b_x2geom[i][j + 1] != 0) ? c_x2geom[i][j + 1] * u[i][j] : 0));

		}

		if (rnorm == 0) { //still zero, we don't need to iterate
			for (i = 0; i < ZMAX; i++) for (j = 0; j < RMAX; j++) u_in[i][j] = 0;
			return 0;
		}
	}
	rnorm = sqrt(rnorm);
	res = sqrt(rsum) / rnorm;
#ifdef DADI_DEBUG
	printf("dadi: res = %g\n", res);
	printf("dadi: rnorm= %g\n", rnorm);
#endif

	if (res < tol_test) return 0;  // no need to iterate

	/*************************************************/
	if (del_t == 0.0) del_t = del_t0; else del_t /= 4;
	del_td = 2.0 * del_t;
	ndiscard = 0;

	/********************/
	/* Begin iteration. */

	for (iter = 0; iter < itermax; iter++)
	{
		/*************************************************/
		/* Copy u into the work array and storage array. */

		for (i = 0; i < ZMAX; i++)
			for (j = 0; j < RMAX; j++) uwork[i][j] = ustor[i][j] = u[i][j];

		/************************************/
		/* Two advances of u via ADI at del_t. */

		adi(u, s, del_t);
		adi(u, s, del_t);

		/*****************************************/
		/* One advance of uwork via ADI at 2*del_t. */

		adi(uwork, s, del_td);

		/*******************************************************/
		/* Calculate test parameter and normalized error.
		   For Dirichlet BCs, no need to worry about boundary
		   points since u,uwork, and ustor should be the same. */

		tptop = tpbot = rsum = 0.0;

		for (i = 1; i < ZMAX; i++)
			for (j = 1; j < RMAX; j++)
			{
				/* Test paramter sums. */
				tptop += (u[i][j] - uwork[i][j]) * (u[i][j] - uwork[i][j]);
				tpbot += (u[i][j] - ustor[i][j]) * (u[i][j] - ustor[i][j]);

				/* Residual terms. */

				dxdxutrm = a_x1geom[i][j] * u[i - 1][j] + b_x1geom[i][j] * u[i][j] + c_x1geom[i][j] * u[i + 1][j];
				dydyutrm = a_x2geom[i][j] * u[i][j - 1] + b_x2geom[i][j] * u[i][j] + c_x2geom[i][j] * u[i][j + 1];

				/* Use the residual as the absolute error and if it is bigger
				   than the stored maximum absolute error, record new maximum
				   absolute error and location.  */
				   /* only include points which are not in structures. */
				errchk = ((b_x2geom[i][j] == 0) ? 0 : dxdxutrm + dydyutrm - s[i][j]);

				/* Residual sums. */
				rsum += errchk * errchk;
			}

		/* Calculate normalized residual. */
		res = sqrt(rsum) / rnorm;
#ifdef DADI_DEBUG
		//printf("dadi: iter= %d, res = %lg del_t=%le\n", iter, res,del_t);
		printf("dadi: iter= %d, res = %g del_t=%e\n", iter, res, del_t);
#endif // DADI_DEBUG
		/* If the residual is less than the tolerance, SUCCESS! */
		if ((res < tol_test) && (iter))
		{
#ifdef DADI_DEBUG
			printf("dadi: iter=%d\n", iter);
#endif// DADI_DEBUG
			for (i = 0; i < ZMAX; i++)
				for (j = 0; j < RMAX; j++)
					u_in[i][j] = (double)u[i][j];

			return(0);
		}

		/* Determine ratio used to find the time step change.  If tpbot
		   is zero but tptop is finite, consider this a case of a large
		   ratio and act accordingly.  DWH does about the same thing
		   except he does NOT discard the solution if tpbot=0. */

		if (tpbot > 0.0) ratio = tptop / tpbot;
		if (tpbot == 0.0) ratio = 1.0;
#ifndef NO_STEP_ADJUST    
		/* Get next time step. */
		if (ratio < 0.02) del_t *= 8.0;
		if (ratio >= 0.02 && ratio < 0.05) del_t *= 4.0;
		if (ratio >= 0.05 && ratio < 0.10) del_t *= 2.0;
		if (ratio >= 0.10 && ratio < 0.30) del_t *= 0.80;
		if (ratio >= 0.30 && ratio < 0.40) del_t *= 0.50;
		if (ratio >= 0.40 && ratio < 0.60) del_t *= 0.25;
		for (i = 0; i < ZMAX; i++)
			for (j = 0; j < RMAX; j++)
				u_in[i][j] = (double)u[i][j];
#endif   

		/* Ratio is too large. */
		if (ratio >= 0.60)
		{
			ndiscard++;
			iter--;
#ifdef DADI_DEBUG
			//  printf("ndiscard= %d, iter=%d step=%lf\n", ndiscard, iter,del_t); 
			printf("ndiscard= %d, iter=%d step=%f\n", ndiscard, iter, del_t);
#endif //DADI_DEBUG

			/* Check if too many discards. */
			if (ndiscard > 20)
			{
				for (i = 0; i < ZMAX; i++)
					for (j = 0; j < RMAX; j++)
						u_in[i][j] = (double)u[i][j];
				del_t = del_t0;
				//		  if(solve(u_in,s,itermax,tol_test))
				printf("Poisson solve FAILURE: dadi: iter= %d, ndiscard>20\n", iter);
				return 1;
			}
			/* Discard by replacing u with what we started with. */
			for (i = 0; i < ZMAX; i++)
				for (j = 0; j < RMAX; j++) u[i][j] = ustor[i][j];

			/* Reduce del_t. */
			del_t /= 8.0;
			//del_t = del_t0;
		}
		del_td = 2 * del_t;
	}
	/* Fail if used up the maximum iterations. */

	printf("Poisson solve FAILURE: dadi:  iter>= %d\n", itermax);

	for (i = 0; i < ZMAX; i++)
		for (j = 0; j < RMAX; j++)
			u_in[i][j] = (double)u[i][j];

	return(2);
}

/***********************************************************************
  Tridiagonal field solver:

  | b0 c0 0                      | | u0 |     | r0 |
  | a1 b1 c1                     | | u1 |     | r1 |
  |      ............            | | .  |  =  | .  |
  |               an-2 bn-2 cn-2 | |un-2|     |rn-2|
  |               0    an-1 bn-1 | |un  |     |rn-1|

  **********************************************************************/

void tridag(double* a, double* b, double* c, double* r, double* utri, double* gam, int n)
{
	register int i;
	double bet;

	/*******************************************/
	/* Decomposition and forward substitution. */

	bet = b[0];
	utri[0] = r[0] / bet;

	for (i = 1; i < n; i++)
	{
		gam[i] = c[i - 1] / bet;
		bet = b[i] - a[i] * gam[i];
		utri[i] = (r[i] - a[i] * utri[i - 1]) / bet;
	}

	/**********************/
	/* Back substitution. */
	for (i = n - 2; i >= 0; i--) utri[i] -= gam[i + 1] * utri[i + 1];
}


//void electric_field()
//{
//	for (int i = 0; i < ZMAX; i++)
//	{
//		for (int j = 0; j < RMAX; j++)
//		{
//			if (btype[i][j] == 1)
//			{
//				Ez[i][j] = -(phi[i + 1][j] - phi[i - 1][j]) / (2 * dz);
//				Er[i][j] = -(phi[i][j + 1] - phi[i][j - 1]) / (2 * dr);
//			}
//			else if (btype[i][j] == LEFT)
//			{
//				Ez[i][j] = (-3 * phi[i][j] + 4 * phi[i + 1][j] - phi[i + 2][j]) / (2 * dz);
//				Er[i][j] = -(phi[i][j + 1] - phi[i][j - 1]) / (2 * dr);
//			}
//			else if (btype[i][j] == RIGHT)
//			{
//				Ez[i][j] = -(-3 * phi[i][j] + 4 * phi[i - 1][j] - phi[i - 2][j]) / (2 * dz);
//				Er[i][j] = -(phi[i][j + 1] - phi[i][j - 1]) / (2 * dr);
//			}
//			else if (btype[i][j] == DOWN)
//			{
//				Ez[i][j] = -(phi[i + 1][j] - phi[i - 1][j]) / (2 * dz);
//				Er[i][j] = (-phi[i][j + 2] + 4 * phi[i][j + 1] - 3 * phi[i][j]) / (2 * dr);
//			}
//			else if (btype[i][j] == UP)
//			{
//				Ez[i][j] = -(phi[i + 1][j] - phi[i - 1][j]) / (2 * dz);
//				Er[i][j] = (phi[i][j - 2] - 4 * phi[i][j - 1] + 3 * phi[i][j]) / (2 * dr);
//			}
//			else if (btype[i][j] == (LEFT + UP))
//			{
//				Ez[i][j] = (-3 * phi[i][j] + 4 * phi[i + 1][j] - phi[i + 2][j]) / (2 * dz);
//				Er[i][j] = (phi[i][j - 2] - 4 * phi[i][j - 1] + 3 * phi[i][j]) / (2 * dr);
//			}
//			else if (btype[i][j] == (LEFT + DOWN))
//			{
//				Ez[i][j] = (-3 * phi[i][j] + 4 * phi[i + 1][j] - phi[i + 2][j]) / (2 * dz);
//				Er[i][j] = (-phi[i][j + 2] + 4 * phi[i][j + 1] - 3 * phi[i][j]) / (2 * dr);
//			}
//			else if (btype[i][j] == (RIGHT + DOWN))
//			{
//				Ez[i][j] = -(-3 * phi[i][j] + 4 * phi[i - 1][j] - phi[i - 2][j]) / (2 * dz);
//				Er[i][j] = (-phi[i][j + 2] + 4 * phi[i][j + 1] - 3 * phi[i][j]) / (2 * dr);
//			}
//			else if (btype[i][j] == (RIGHT + UP))
//			{
//				Ez[i][j] = -(-3 * phi[i][j] + 4 * phi[i - 1][j] - phi[i - 2][j]) / (2 * dz);
//				Er[i][j] = (phi[i][j - 2] - 4 * phi[i][j - 1] + 3 * phi[i][j]) / (2 * dr);
//
//			}
//			else if (btype[i][j] == 0)
//			{
//				Ez[i][j] = 0;
//				Er[i][j] = 0;
//
//			}
//
//		}
//	}
//}

void electric_field()
{
	int i, j, k;
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			if (btype[i][j] == 0 || btype[i][j] == 110) continue;

			if (i > 0 &&j > 0 && i< (ZMAX - 1) &&j<(RMAX - 1))
			{
				Ez[i][j] = -(phi[i + 1][j] - phi[i - 1][j]) / (2 * dz);
				Er[i][j] = -(phi[i][j + 1] - phi[i][j - 1]) / (2 * dr);
			}
			else if (i == 0 && j > 0  && j < (RMAX - 1))
			{
				Ez[i][j] = -(-3 * phi[i][j] + 4 * phi[i + 1][j] - phi[i + 2][j]) / (2 * dz);
				Er[i][j] = -(phi[i][j + 1] - phi[i][j - 1]) / (2 * dr);
			}
			else if (i == (ZMAX - 1) && j > 0 && j < (RMAX - 1))
			{
				Ez[i][j] = -(3 * phi[i][j] - 4 * phi[i - 1][j] + phi[i - 2][j]) / (2 * dz);
				Er[i][j] = -(phi[i][j + 1] - phi[i][j - 1]) / (2 * dr);
			}
			else if (i > 0 && j == 0 && i < (ZMAX - 1))
			{
				Ez[i][j] = -(phi[i + 1][j] - phi[i - 1][j]) / (2 * dz);
				Er[i][j] = -(-phi[i][j + 2] + 4 * phi[i][j + 1] - 3 * phi[i][j]) / (2 * dr);
			}
			else if (i > 0 && j == (RMAX - 1) && i < (ZMAX - 1))
			{
				Ez[i][j] = -(phi[i + 1][j] - phi[i - 1][j]) / (2 * dz);
				Er[i][j] = -(phi[i][j - 2] - 4 * phi[i][j - 1] + 3 * phi[i][j]) / (2 * dr);
			}
			else if (i == 0 && j == (RMAX - 1))
			{
				Ez[i][j] = -(-3 * phi[i][j] + 4 * phi[i + 1][j] - phi[i + 2][j]) / (2 * dz);
				Er[i][j] = -(phi[i][j - 2] - 4 * phi[i][j - 1] + 3 * phi[i][j]) / (2 * dr);
			}
			else if (i == 0 && j == 0)
			{
				Ez[i][j] = -(-3 * phi[i][j] + 4 * phi[i + 1][j] - phi[i + 2][j]) / (2 * dz);
				Er[i][j] = -(-3 * phi[i][j] + 4 * phi[i][j + 1] - phi[i][j + 2]) / (2 * dr);
			}
			else if (i == (ZMAX - 1) && j == 0)
			{
				Ez[i][j] = -(3 * phi[i][j] - 4 * phi[i - 1][j] + phi[i - 2][j]) / (2 * dz);
				Er[i][j] = -(-3 * phi[i][j] + 4 * phi[i][j + 1] - phi[i][j + 2]) / (2 * dr);
			}
			else if (i == (ZMAX - 1) && j == (RMAX - 1))
			{
				Ez[i][j] = -(3 * phi[i][j] - 4 * phi[i - 1][j] + phi[i - 2][j]) / (2 * dz);
				Er[i][j] = -(3 * phi[i][j] - 4 * phi[i][j - 1] + phi[i][j - 2]) / (2 * dr);

			}
		}
	}

	//_for(k, 0, BND_NUM)
	//{
	//	if (boundary_array[k].physics_type == DIELECTRIC_SURFACE_BOUNDARY)
	//	{
	//		double dr = boundary_array[k].end.r - boundary_array[k].start.r;
	//		double dz = boundary_array[k].end.z - boundary_array[k].start.z;

	//		if (dr > dz)
	//		{
	//			double ins_z = dz / dr;
	//			i = 0;
	//			_feq(j, boundary_array[k].start.r, boundary_array[k].end.r)
	//			{
	//				int ti = (int)(boundary_array[k].start.z + ceil(i * ins_z));
	//				if (btype[ti][j] == UP || btype[ti][j] == DOWN)
	//				{
	//					Er[ti][j] = -(MPDT[ti][j].ni - MPDT[ti][j].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				else if (btype[ti][j] == LEFT || btype[ti][j] == RIGHT)
	//				{
	//					Ez[ti][j] = -(MPDT[ti][j].ni - MPDT[ti][j].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				else if (btype[ti][j] == (LEFT + UP))
	//				{
	//					Er[ti][j] = -0.7071 * (MPDT[ti][j].ni - MPDT[ti][j].ne) * QE / EPS_0 * EPS_PLA;
	//					Ez[ti][j] = -0.7071 * (MPDT[ti][j].ni - MPDT[ti][j].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				else if (btype[ti][j] == (LEFT + DOWN))
	//				{
	//					Er[ti][j] = -0.7071 * (MPDT[ti][j].ni - MPDT[ti][j].ne) * QE / EPS_0 * EPS_PLA;
	//					Ez[ti][j] = -0.7071 * (MPDT[ti][j].ni - MPDT[ti][j].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				else if (btype[ti][j] == (RIGHT + DOWN))
	//				{
	//					Er[ti][j] = -0.7071 * (MPDT[ti][j].ni - MPDT[ti][j].ne) * QE / EPS_0 * EPS_PLA;
	//					Ez[ti][j] = -0.7071 * (MPDT[ti][j].ni - MPDT[ti][j].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				else if (btype[ti][j] == (RIGHT + UP))
	//				{
	//					Er[ti][j] = -0.7071 * (MPDT[ti][j].ni - MPDT[ti][j].ne) * QE / EPS_0 * EPS_PLA;
	//					Ez[ti][j] = -0.7071 * (MPDT[ti][j].ni - MPDT[ti][j].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				i++;
	//			}
	//		}
	//		else
	//		{
	//			double ins_r = dr / dz;
	//			j = 0;
	//			_feq(i, boundary_array[k].start.z, boundary_array[k].end.z)
	//			{
	//				int tj = (int)(boundary_array[k].start.r + ceil(j * ins_r));
	//				if (btype[i][tj] == UP || btype[i][tj] == DOWN)
	//				{
	//					Er[i][tj] = -(MPDT[i][tj].ni - MPDT[i][tj].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				else if (btype[i][tj] == LEFT || btype[i][tj] == RIGHT)
	//				{
	//					Ez[i][tj] = -(MPDT[i][tj].ni - MPDT[i][tj].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				else if (btype[i][tj] == (LEFT + UP))
	//				{
	//					Er[i][tj] = -0.7071 * (MPDT[i][tj].ni - MPDT[i][tj].ne) * QE / EPS_0 * EPS_PLA;
	//					Ez[i][tj] = -0.7071 * (MPDT[i][tj].ni - MPDT[i][tj].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				else if (btype[i][tj] == (LEFT + DOWN))
	//				{
	//					Er[i][tj] = -0.7071 * (MPDT[i][tj].ni - MPDT[i][tj].ne) * QE / EPS_0 * EPS_PLA;
	//					Ez[i][tj] = -0.7071 * (MPDT[i][tj].ni - MPDT[i][tj].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				else if (btype[i][tj] == (RIGHT + DOWN))
	//				{
	//					Er[i][tj] = -0.7071 * (MPDT[i][tj].ni - MPDT[i][tj].ne) * QE / EPS_0 * EPS_PLA;
	//					Ez[i][tj] = -0.7071 * (MPDT[i][tj].ni - MPDT[i][tj].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				else if (btype[i][tj] == (RIGHT + UP))
	//				{
	//					Er[i][tj] = -0.7071 * (MPDT[i][tj].ni - MPDT[i][tj].ne) * QE / EPS_0 * EPS_PLA;
	//					Ez[i][tj] = -0.7071 * (MPDT[i][tj].ni - MPDT[i][tj].ne) * QE / EPS_0 * EPS_PLA;
	//				}
	//				j++;
	//			}
	//		}
	//	}
	//}
}

void potential_boundary()
{
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			if (btype[i][j] == 1)
			{

			}
			else if (btype[i][j] == LEFT)
			{
				phi[i][j] = phi[i + 1][j];
			}
			else if (btype[i][j] == RIGHT)
			{
				phi[i][j] = phi[i - 1][j];
			}
			else if (btype[i][j] == UP)
			{
				phi[i][j] = phi[i][j - 1];
			}
			else if (btype[i][j] == DOWN)
			{
				phi[i][j] = phi[i][j + 1];
			}
			else if (btype[i][j] == (LEFT + UP))
			{
				phi[i][j] = phi[i + 1][j - 1];
			}
			else if (btype[i][j] == (LEFT + DOWN))
			{
				phi[i][j] = phi[i + 1][j + 1];
			}
			else if (btype[i][j] == (RIGHT + DOWN))
			{
				phi[i][j] = phi[i - 1][j + 1];
			}
			else if (btype[i][j] == (RIGHT + UP))
			{
				phi[i][j] = phi[i - 1][j - 1];
			}
			else if (btype[i][j] == 0)
			{

			}
		}
	}
}

void potential_solve()
{
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			if (MPDT[i][j].peq > MPDT[i][j].neq)
			{
				MPDT[i][j].peq -= MPDT[i][j].neq;
				MPDT[i][j].neq = 0;
			}
			else
			{
				MPDT[i][j].neq -= MPDT[i][j].peq;
				MPDT[i][j].peq = 0;
			}

			//rho[i][j] = -(MPDT[i][j].ni - MPDT[i][j].ne) * QE / EPS_0 * EPS_PLA;
			rho[i][j] = -(MPDT[i][j].peq - MPDT[i][j].neq) * QE / EPS_0 * EPS_PLA;
			rou[i][j] = (MPDT[i][j].peq - MPDT[i][j].neq);
		}
	}
	solve(phi, rho, 50, 0.01);

	max_phi = 0;
	for (int i = 0; i < ZMAX; i++)
	{
		for (int j = 0; j < RMAX; j++)
		{
			if (abs(phi[i][j]) > max_phi)
			{
				max_phi = abs(phi[i][j]);
			}
		}
	}
	printf("max phi = %lf\n", max_phi);
	//solveGS();
	//potential_boundary();

	electric_field();
}

void mag_phi()
{
	register int i, j;
	for (i = 1; i < nz - 1; i++)
	{
		for (j = 1; j < nr - 1; j++)
		{
			

			if (i > 0 && j > 0 && i < (ZMAX - 1) && j < (RMAX - 1))
			{
				btheta[i][j] += dt * ((Ez[i][j + 1] - Ez[i][j - 1]) / (2 * dr) - (Er[i + 1][j] - Er[i - 1][j]) / (2 * dz));
			}
			else if (i == 0 && j > 0 && j < (RMAX - 1))
			{
				btheta[i][j] += dt * ((Ez[i][j + 1] - Ez[i][j - 1]) / (2 * dr) - (-3 * Er[i][j] + 4 * Er[i + 1][j] - Er[i + 2][j]) / (2 * dz));
			}
			else if (i == (ZMAX - 1) && j > 0 && j < (RMAX - 1))
			{
				btheta[i][j] += dt * ((Ez[i][j + 1] - Ez[i][j - 1]) / (2 * dr) - (3 * Er[i][j] - 4 * Er[i + 1][j] + Er[i - 2][j]) / (2 * dz));
			}
			else if (i > 0 && j == 0 && i < (ZMAX - 1))
			{
				btheta[i][j] += dt * ((-3 * Ez[i][j] + 4 * Ez[i][j + 1] - Ez[i][j + 2]) / (2 * dr) - (Er[i + 1][j] - Er[i - 1][j]) / (2 * dz));
			}
			else if (i > 0 && j == (RMAX - 1) && i < (ZMAX - 1))
			{
				btheta[i][j] += dt * (( 3 * Ez[i][j] - 4 * Ez[i][j - 1] + Ez[i][j - 2]) / (2 * dr) - (Er[i + 1][j] - Er[i - 1][j]) / (2 * dz));
			}
			else if (i == 0 && j == (RMAX - 1))
			{
				btheta[i][j] += dt * ((3 * Ez[i][j] - 4 * Ez[i][j - 1] + Ez[i][j - 2]) / (2 * dr) - (-3 * Er[i][j] + 4 * Er[i + 1][j] - Er[i + 2][j]) / (2 * dz));
			}
			else if (i == 0 && j == 0)
			{
				btheta[i][j] += dt * ((-3 * Ez[i][j] + 4 * Ez[i][j + 1] - Ez[i][j + 2]) / (2 * dr) - (-3 * Er[i][j] + 4 * Er[i + 1][j] - Er[i + 2][j]) / (2 * dz));
			}
			else if (i == (ZMAX - 1) && j == 0)
			{
				btheta[i][j] += dt * ((-3 * Ez[i][j] + 4 * Ez[i][j + 1] - Ez[i][j + 2]) / (2 * dr) - (3 * Er[i][j] - 4 * Er[i + 1][j] + Er[i - 2][j]) / (2 * dz));
			}
			else if (i == (ZMAX - 1) && j == (RMAX - 1))
			{
				btheta[i][j] += dt * ((3 * Ez[i][j] - 4 * Ez[i][j - 1] + Ez[i][j - 2]) / (2 * dr) - (-3 * Er[i][j] + 4 * Er[i + 1][j] - Er[i + 2][j]) / (2 * dz));

			}
		}
	}


}