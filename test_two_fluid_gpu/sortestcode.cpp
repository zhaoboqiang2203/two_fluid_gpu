#include "MHD.h"
using namespace std;

void outputcsv()
{
	ofstream out("outfile.csv");

	for (int i = 0; i < nz; i++)
	{
		for (int j = 0; j < nr; j++)
		{
			out << phi[i][j] << ",";
		}
		out << endl;
	}

}

void test_sor_code()
{
	for (int i = 0; i < nz; i++)
	{
		for (int j = 0; j < nr; j++)
		{
			phi[i][j] = 0;

		}
	}

	for (int i = 40; i < 60; i++)
	{
		for (int j = 40; j < 60; j++)
		{
			rho[i][j] = 1e-17;
		}
	}

	solveGS();

	outputcsv();
}