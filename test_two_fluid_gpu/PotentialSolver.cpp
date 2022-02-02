#include "MHD.h" 

int solveGS()
{
    //references to avoid having to write world.phi1
    //rou contains only ion contribution

    //TODO 根据实际网格大小修改
    double idx2 = dz;
    double idy2 = dr;
    double phi0 = 0;		//reference plasma potential
    double n0 = 0;		//reference electron density
    double Te0 = 1.5;		//reference electron temperature in eV
    double L2 = 0;			//norm
    double tolerance = 10e-8;
    int converged = 0;
    int max_solver_it = 1000;
    /*solve potential*/
    for (unsigned it = 0; it < max_solver_it; it++)
    {
        //printf(" it = %d \n", it);
        for (int i = 0; i < nz; i++)
            for (int j = 0; j < nr; j++)
            {
                /*skip over solid (fixed) nodes = Dirichlet boundaries*/
                //if (world.object_id[i][j]>0) continue;

                //TODO 根据world的标记修改
                //if (world[i][j] == 0)
                //	phi1[i][j] = phi1[i+1][j];
                //else if (i==world.ni-1)
                //	phi1[i][j] = phi1[i-1][j];
                //else if (j==0)
                //	phi1[i][j] = phi1[i][j+1];
                //else if (j==world.nj-1)
                //	phi1[i][j] = phi1[i][j-1];
                //else if (k==0)
                //	phi1[i][j] = phi1[i][j][k+1];
                //else if (k==world.nk-1)
                //	phi1[i][j] = phi1[i][j][k-1];
                //else {	//standard internal open node

                    //evaluate electron density from the Boltzmann relationshp


                double tne = n0 * exp((phi1[i][j] - phi0) / Te0);

                double phi_new = ((rou[i][j] - QE * tne) / EPS_0 +
                    idx2 * (phi1[i - 1][j] + phi1[i + 1][j]) +
                    idy2 * (phi1[i][j - 1] + phi1[i][j + 1])) / (2 * idx2 + 2 * idy2);

                /*SOR*/
                phi1[i][j] = phi1[i][j] + 1.4 * (phi_new - phi1[i][j]);
            }

        /*check for convergence*/
        if (it % 25 == 0)
        {
            double sum = 0;
            for (int i = 0; i < nz; i++)
                for (int j = 0; j < nr; j++)

                {
                    /*skip over solid (fixed) nodes*/
                    //if (world.object_id[i][j]>0) continue;

                    double R = 0;
                    //if (i==0)
                    //	R = phi1[i][j] - phi1[i+1][j];
                    //else if (i==ni-1)
                    //	R = phi1[i][j] - phi1[i-1][j];
                    //else if (j==0)
                    //	R = phi1[i][j] - phi1[i][j+1];
                    //else if (j==world.nj-1)
                    //	R = phi1[i][j] - phi1[i][j-1];
                    //else {
                            //evaluate electron density from the Boltzmann relationshp
                    double tne = n0 * exp((phi1[i][j] - phi0) / Te0);
                    R = -phi1[i][j] * (2 * idx2 + 2 * idy2) +
                        (rou[i][j] - QE * tne) / EPS_0 +
                        idx2 * (phi1[i - 1][j] + phi1[i + 1][j]) +
                        idy2 * (phi1[i][j - 1] + phi1[i][j + 1]);

                    sum += R * R;
                }

            L2 = sqrt(sum / (nr * nz));
            printf("L2 = %lf \n", L2);
            if (L2 < tolerance) { converged = true; break; }
        }
    }

    if (!converged) printf("GS failed to converge, L2 = %lf \n", L2);
    return converged;
}