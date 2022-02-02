#include "testcode.h"

void test_potential()
{
	double res[ZMAX][RMAX];
	double test_s[ZMAX][RMAX];
	//f = cos()

	for(int i=0;i<RMAX;i++)
		for (int j = 0; j < RMAX; j++)
		{
			res[i][j] = 0;
			test_s[i][j] = 0;
		}

	for (int i = 60; i < 70; i++)
		for (int j = 60; j < 70; j++)
		{
			test_s[i][j] = 1;
		}
	//f = cos();
	matrix_to_csv((double**)test_s, ZMAX, RMAX, RMAX, (char*)(".\\output\\test_s.csv"));
	solve(res, test_s, 1000, 0.1);
	matrix_to_csv((double**)res, ZMAX, RMAX, RMAX, (char*)(".\\output\\testphi.csv"));
}
