#include "fluid.h"
using namespace std;

float Une[ZMAX][RMAX];
float Uni[ZMAX][RMAX];
float Uver[ZMAX][RMAX];
float Uvetheta[ZMAX][RMAX];
float Uvez[ZMAX][RMAX];
float Uvir[ZMAX][RMAX];
float Uvitheta[ZMAX][RMAX];
float Uviz[ZMAX][RMAX];
float Br[ZMAX][RMAX];
float Btheta[ZMAX][RMAX];
float Bz[ZMAX][RMAX];
float Pe[ZMAX][RMAX];
float Pi[ZMAX][RMAX];

float Une_bar[ZMAX][RMAX];
float Uni_bar[ZMAX][RMAX];
float Uver_bar[ZMAX][RMAX];
float Uvetheta_bar[ZMAX][RMAX];
float Uvez_bar[ZMAX][RMAX];
float Uvir_bar[ZMAX][RMAX];
float Uvitheta_bar[ZMAX][RMAX];
float Uviz_bar[ZMAX][RMAX];
float Br_bar[ZMAX][RMAX];
float Btheta_bar[ZMAX][RMAX];
float Bz_bar[ZMAX][RMAX];
float Pe_bar[ZMAX][RMAX];
float Pi_bar[ZMAX][RMAX];

float F1[ZMAX][RMAX];
float F2[ZMAX][RMAX];
float F3[ZMAX][RMAX];
float F4[ZMAX][RMAX];
float F5[ZMAX][RMAX];
float F6[ZMAX][RMAX];
float F7[ZMAX][RMAX];
float F8[ZMAX][RMAX]; 
float F9[ZMAX][RMAX];
float F10[ZMAX][RMAX];
float F11[ZMAX][RMAX];
float F12[ZMAX][RMAX]; 
float F13[ZMAX][RMAX];


//__global__ void Kernel_test

int row = 386, col = 1;

//cudaError_t Maccormack_(int* c, const int* a, const int* b, unsigned int size)
//{
//	int* dev_a = 0;
//	int* dev_b = 0;
//	int* dev_c = 0;
//	cudaError_t cudaStatus;
//
//	// Choose which GPU to run on, change this on a multi-GPU system.
//	cudaStatus = cudaSetDevice(0);
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
//		goto Error;
//	}
//
//	// Allocate GPU buffers for three vectors (two input, one output)    .
//	cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMalloc failed!");
//		goto Error;
//	}
//
//	cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMalloc failed!");
//		goto Error;
//	}
//
//	cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMalloc failed!");
//		goto Error;
//	}
//
//	// Copy input vectors from host memory to GPU buffers.
//	cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMemcpy failed!");
//		goto Error;
//	}
//
//	cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMemcpy failed!");
//		goto Error;
//	}
//
//	// Launch a kernel on the GPU with one thread for each element.
//	addKernel << <1, size >> > (dev_c, dev_a, dev_b);
//
//	// Check for any errors launching the kernel
//	cudaStatus = cudaGetLastError();
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
//		goto Error;
//	}
//
//	// cudaDeviceSynchronize waits for the kernel to finish, and returns
//	// any errors encountered during the launch.
//	cudaStatus = cudaDeviceSynchronize();
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
//		goto Error;
//	}
//
//	// Copy output vector from GPU buffer to host memory.
//	cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
//	if (cudaStatus != cudaSuccess) {
//		fprintf(stderr, "cudaMemcpy failed!");
//		goto Error;
//	}
//
//Error:
//	cudaFree(dev_c);
//	cudaFree(dev_a);
//	cudaFree(dev_b);
//
//	return cudaStatus;
//}


void electron_flow_v2()
{
	struct _U Qr;
	struct _U Qz;
	register int i, j, k;
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{

			U[i][j] = cal_u(i, j);
			Fr[i][j] = cal_fr(U[i][j]);
			Fz[i][j] = cal_fz(U[i][j]);
			s[i][j] = cal_s(U[i][j]);

			U_bar[i][j].z = i;
			U_bar[i][j].r = j;
			U_bar2[i][j].z = i;
			U_bar2[i][j].r = j;

		}
	}

	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				for (k = 0; k < 13; k++)
				{
					U_bar[i][j].u[k] = U[i][j].u[k] - dt / dr * (Fr[i][j + 1].f[k] - Fr[i][j].f[k]) - dt / dz * (Fz[i + 1][j].f[k] - Fz[i][j].f[k]) + dt * s[i][j].f[k];
				}
#ifdef FLUID_DEBUG
				if (i == row && j == col)
				{
					printf("U[i][j].u1 = %.5e\n", U[i][j].u[0]);
					printf("Fr[i][j + 1].f1 = %.5e\n", Fr[i][j + 1].f[0]);
					printf("Fr[i][j].f1 = %.5e\n", Fr[i][j].f[0]);
					printf("Fz[i + 1][j].f1 = %.5e\n", Fz[i + 1][j].f[0]);
					printf("Fz[i][j].f1 = %.5e\n", Fz[i][j].f[0]);
					printf("s[i][j].f1 = %.5e\n", s[i][j].f[0]);
					printf("U_bar[i][j].u1 = %.5e\n", U_bar[i][j].u[0]);

					printf("U[i][j].u3 = %.5e\n", U[i][j].u[2]);
					printf("Fr[i][j + 1].f3 = %.5e\n", Fr[i][j + 1].f[2]);
					printf("Fr[i][j].f3 = %.5e\n", Fr[i][j].f[2]);
					printf("Fz[i + 1][j].f3 = %.5e\n", Fz[i + 1][j].f[2]);
					printf("Fz[i][j].f3 = %.5e\n", Fz[i][j].f[2]);
					printf("s[i][j].f3 = %.5e\n", s[i][j].f[2]);
					printf("U_bar[i][j].u3 = %.5e\n", U_bar[i][j].u[2]);
				}
#endif// FLUID_DEBUG
			}
			else if (btype[i][j] == 0)
			{
				for (k = 0; k < 13; k++)
				{
					U_bar[i][j].u[k] = 0;
				}
			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				for (k = 0; k < 13; k++)
				{
					//U_bar[i][j].u[k] = U[i][j].u[k] - dt / dr * (Fr[i][j + 1].f[k] - Fr[i][j].f[k]) - dt / dz * (Fz[i + 1][j].f[k] - Fz[i][j].f[k]) + dt * s[i][j].f[k];
					U_bar[i][j].u[k] = U[i][j].u[k] -  dt / dz * (Fz[i + 1][j].f[k] - Fz[i][j].f[k]) + dt * s[i][j].f[k];
				}
			}
			else
			{
				for (k = 0; k < 13; k++)
				{
					U_bar[i][j].u[k] = U[i][j].u[k];
				}
			}

			Fr_bar[i][j] = cal_fr(U_bar[i][j]);
			Fz_bar[i][j] = cal_fz(U_bar[i][j]);
			s_bar[i][j] = cal_s(U_bar[i][j]);
		}
	}

	//矫正步
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				for (k = 0; k < 13; k++)
				{
					U_bar2[i][j].u[k] = 0.5 * (U[i][j].u[k] + U_bar[i][j].u[k] - dt / dr * (Fr_bar[i][j].f[k] - Fr_bar[i][j - 1].f[k]) - dt / dz * (Fz_bar[i][j].f[k] - Fz_bar[i - 1][j].f[k]) + dt * s_bar[i][j].f[k]);
				}
#ifdef FLUID_DEBUG
				if (i == row && j == col)
				{
					printf("Fr_bar[i][j].f[0] = %.5e\n", Fr_bar[i][j].f[0]);
					printf("Fr_bar[i][j - 1].f[0] = %.5e\n", Fr_bar[i][j - 1].f[0]);
					printf("Fz_bar[i][j].f[0] = %.5e\n", Fz_bar[i][j].f[0]);
					printf("Fz_bar[i - 1][j].f[0] = %.5e\n", Fz_bar[i - 1][j].f[0]);
					printf("s_bar[i][j].f[0] = %.5e\n", s_bar[i][j].f[0]);
					printf("U_bar2[i][j].u[0] = %.5e\n", U_bar2[i][j].u[0]);

					printf("Fr_bar[i][j].f3 = %.5e\n", Fr_bar[i][j].f[2]);
					printf("Fr_bar[i][j - 1].f3 = %.5e\n", Fr_bar[i][j - 1].f[2]);
					printf("Fz_bar[i][j].f3 = %.5e\n", Fz_bar[i][j].f[2]);
					printf("Fz_bar[i - 1][j].f3 = %.5e\n", Fz_bar[i - 1][j].f[2]);
					printf("s_bar[i][j].f3 = %.5e\n", s_bar[i][j].f[2]);
					printf("U_bar2[i][j].u3 = %.5e\n", U_bar2[i][j].u[2]);

				}
#endif// FLUID_DEBUG				

			}
			else if (btype[i][j] == 0)
			{
				for (k = 0; k < 13; k++)
				{
					U_bar2[i][j].u[k] = 0;
				}
			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				for (k = 0; k < 13; k++)
				{
					//U_bar2[i][j].u[k] = 0.5 * (U[i][j].u[k] + U_bar[i][j].u[k] - dt / dr * (Fr_bar[i][j + 1].f[k] - Fr_bar[i][j].f[k]) - dt / dz * (Fz_bar[i][j].f[k] - Fz_bar[i - 1][j].f[k]) + dt * s_bar[i][j].f[k]);
					U_bar2[i][j].u[k] = 0.5 * (U[i][j].u[k] + U_bar[i][j].u[k] - dt / dz * (Fz_bar[i][j].f[k] - Fz_bar[i - 1][j].f[k]) + dt * s_bar[i][j].f[k]);
				}
			}
			else
			{
				for (k = 0; k < 13; k++)
				{
					U_bar2[i][j].u[k] = U_bar[i][j].u[k];
				}
			}

		}
	}

	//计算				

	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				Qz = arti_vis(MPDT[i + 1][j], MPDT[i][j], MPDT[i - 1][j]);
				Qr = arti_vis(MPDT[i][j + 1], MPDT[i][j], MPDT[i][j - 1]);

				U[i][j].u[0] = U_bar2[i][j].u[0] + Qr.u[0] / 4 * (MPDT[i][j + 1].ne - 2 * MPDT[i][j].ne + MPDT[i][j - 1].ne) + Qz.u[0] / 4 * (MPDT[i + 1][j].ne - 2 * MPDT[i][j].ne + MPDT[i - 1][j].ne);
				U[i][j].u[1] = U_bar2[i][j].u[1] + Qr.u[1] / 4 * (MPDT[i][j + 1].ni - 2 * MPDT[i][j].ni + MPDT[i][j - 1].ni) + Qz.u[1] / 4 * (MPDT[i + 1][j].ni - 2 * MPDT[i][j].ni + MPDT[i - 1][j].ni);
				U[i][j].u[2] = U_bar2[i][j].u[2] + Qr.u[2] / 4 * (MPDT[i][j + 1].ver - 2 * MPDT[i][j].ver + MPDT[i][j - 1].ver) + Qz.u[2] / 4 * (MPDT[i + 1][j].ver - 2 * MPDT[i][j].ver + MPDT[i - 1][j].ver);
				U[i][j].u[3] = U_bar2[i][j].u[3] + Qr.u[3] / 4 * (MPDT[i][j + 1].vetheta - 2 * MPDT[i][j].vetheta + MPDT[i][j - 1].vetheta) + Qz.u[3] / 4 * (MPDT[i + 1][j].vetheta - 2 * MPDT[i][j].vetheta + MPDT[i - 1][j].vetheta);
				U[i][j].u[4] = U_bar2[i][j].u[4] + Qr.u[4] / 4 * (MPDT[i][j + 1].vez - 2 * MPDT[i][j].vez + MPDT[i][j - 1].vez) + Qz.u[4] / 4 * (MPDT[i + 1][j].vez - 2 * MPDT[i][j].vez + MPDT[i - 1][j].vez);
				U[i][j].u[5] = U_bar2[i][j].u[5] + Qr.u[5] / 4 * (MPDT[i][j + 1].vir - 2 * MPDT[i][j].vir + MPDT[i][j - 1].vir) + Qz.u[5] / 4 * (MPDT[i + 1][j].vir - 2 * MPDT[i][j].vir + MPDT[i - 1][j].vir);
				U[i][j].u[6] = U_bar2[i][j].u[6] + Qr.u[6] / 4 * (MPDT[i][j + 1].vitheta - 2 * MPDT[i][j].vitheta + MPDT[i][j - 1].vitheta) + Qz.u[6] / 4 * (MPDT[i + 1][j].vitheta - 2 * MPDT[i][j].vitheta + MPDT[i - 1][j].vitheta);
				U[i][j].u[7] = U_bar2[i][j].u[7] + Qr.u[7] / 4 * (MPDT[i][j + 1].viz - 2 * MPDT[i][j].viz + MPDT[i][j - 1].viz) + Qz.u[7] / 4 * (MPDT[i + 1][j].viz - 2 * MPDT[i][j].viz + MPDT[i - 1][j].viz);
				U[i][j].u[8] = U_bar2[i][j].u[8] + Qr.u[8] / 4 * (MPDT[i][j + 1].br - 2 * MPDT[i][j].br + MPDT[i][j - 1].br) + Qz.u[8] / 4 * (MPDT[i + 1][j].br - 2 * MPDT[i][j].br + MPDT[i - 1][j].br);
				U[i][j].u[9] = U_bar2[i][j].u[9] + Qr.u[9] / 4 * (MPDT[i][j + 1].btheta - 2 * MPDT[i][j].btheta + MPDT[i][j - 1].btheta) + Qz.u[9] / 4 * (MPDT[i + 1][j].btheta - 2 * MPDT[i][j].btheta + MPDT[i - 1][j].btheta);
				U[i][j].u[10] = U_bar2[i][j].u[10] + Qr.u[10] / 4 * (MPDT[i][j + 1].bz - 2 * MPDT[i][j].bz + MPDT[i][j - 1].bz) + Qz.u[10] / 4 * (MPDT[i + 1][j].bz - 2 * MPDT[i][j].bz + MPDT[i - 1][j].bz);
				U[i][j].u[11] = U_bar2[i][j].u[11] + Qr.u[11] / 4 * (MPDT[i][j + 1].pe - 2 * MPDT[i][j].pe + MPDT[i][j - 1].pe) + Qz.u[11] / 4 * (MPDT[i + 1][j].pe - 2 * MPDT[i][j].pe + MPDT[i - 1][j].pe);
				U[i][j].u[12] = U_bar2[i][j].u[12] + Qr.u[12] / 4 * (MPDT[i][j + 1].pi - 2 * MPDT[i][j].pi + MPDT[i][j - 1].pi) + Qz.u[12] / 4 * (MPDT[i + 1][j].pi - 2 * MPDT[i][j].pi + MPDT[i - 1][j].pi);

#ifdef FLUID_DEBUG
				if (i == row && j == col)
				{
					//cout << "final U[i][j].u1 = " << U[i][j].u1 << endl;
					printf("final U[i][j].u1 = %.5e\n", U[i][j].u[0]);
					printf("U_bar2[i][j + 1].u1 = %.5e\n", U_bar2[i][j + 1].u[0]);
					printf("U_bar2[i][j - 1].u1 = %.5e\n", U_bar2[i][j - 1].u[0]);
					printf("U_bar2[i + 1][j].u1 = %.5e\n", U_bar2[i + 1][j].u[0]);
					printf("U_bar2[i - 1][j].u1 = %.5e\n", U_bar2[i - 1][j].u[0]);

					printf("final U[i][j].u3 = %.5e\n", U[i][j].u[2]);
					printf("U_bar2[i][j + 1].u3 = %.5e\n", U_bar2[i][j + 1].u[2]);
					printf("U_bar2[i][j - 1].u3 = %.5e\n", U_bar2[i][j - 1].u[2]);
					printf("U_bar2[i + 1][j].u3 = %.5e\n", U_bar2[i + 1][j].u[2]);
					printf("U_bar2[i - 1][j].u3 = %.5e\n", U_bar2[i - 1][j].u[2]);

				}
#endif// FLUID_DEBUG	

			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				for (k = 0; k < 13; k++)
				{
					Qz = arti_vis(MPDT[i + 1][j], MPDT[i][j], MPDT[i - 1][j]);
					Qr = arti_vis(MPDT[i][j + 2], MPDT[i][j + 1], MPDT[i][j]);

					U[i][j].u[0] = U_bar2[i][j].u[0] + Qr.u[0] / 2 * (MPDT[i][j + 1].ne - MPDT[i][j].ne) + Qz.u[0] / 4 * (MPDT[i + 1][j].ne - 2 * MPDT[i][j].ne + MPDT[i - 1][j].ne);
					U[i][j].u[1] = U_bar2[i][j].u[1] + Qr.u[1] / 2 * (MPDT[i][j + 1].ni - MPDT[i][j].ni) + Qz.u[1] / 4 * (MPDT[i + 1][j].ni - 2 * MPDT[i][j].ni + MPDT[i - 1][j].ni);
					U[i][j].u[2] = U_bar2[i][j].u[2] + Qr.u[2] / 2 * (MPDT[i][j + 1].ver - MPDT[i][j].ver) + Qz.u[2] / 4 * (MPDT[i + 1][j].ver - 2 * MPDT[i][j].ver + MPDT[i - 1][j].ver);
					U[i][j].u[3] = U_bar2[i][j].u[3] + Qr.u[3] / 2 * (MPDT[i][j + 1].vetheta - MPDT[i][j].vetheta) + Qz.u[3] / 4 * (MPDT[i + 1][j].vetheta - 2 * MPDT[i][j].vetheta + MPDT[i - 1][j].vetheta);
					U[i][j].u[4] = U_bar2[i][j].u[4] + Qr.u[4] / 2 * (MPDT[i][j + 1].vez - MPDT[i][j].vez) + Qz.u[4] / 4 * (MPDT[i + 1][j].vez - 2 * MPDT[i][j].vez + MPDT[i - 1][j].vez);
					U[i][j].u[5] = U_bar2[i][j].u[5] + Qr.u[5] / 2 * (MPDT[i][j + 1].vir - MPDT[i][j].vir) + Qz.u[5] / 4 * (MPDT[i + 1][j].vir - 2 * MPDT[i][j].vir + MPDT[i - 1][j].vir);
					U[i][j].u[6] = U_bar2[i][j].u[6] + Qr.u[6] / 2 * (MPDT[i][j + 1].vitheta - MPDT[i][j].vitheta) + Qz.u[6] / 4 * (MPDT[i + 1][j].vitheta - 2 * MPDT[i][j].vitheta + MPDT[i - 1][j].vitheta);
					U[i][j].u[7] = U_bar2[i][j].u[7] + Qr.u[7] / 2 * (MPDT[i][j + 1].viz - MPDT[i][j].viz) + Qz.u[7] / 4 * (MPDT[i + 1][j].viz - 2 * MPDT[i][j].viz + MPDT[i - 1][j].viz);
					U[i][j].u[8] = U_bar2[i][j].u[8] + Qr.u[8] / 2 * (MPDT[i][j + 1].br - MPDT[i][j].br) + Qz.u[8] / 4 * (MPDT[i + 1][j].br - 2 * MPDT[i][j].br + MPDT[i - 1][j].br);
					U[i][j].u[9] = U_bar2[i][j].u[9] + Qr.u[9] / 2 * (MPDT[i][j + 1].btheta - MPDT[i][j].btheta) + Qz.u[9] / 4 * (MPDT[i + 1][j].btheta - 2 * MPDT[i][j].btheta + MPDT[i - 1][j].btheta);
					U[i][j].u[10] = U_bar2[i][j].u[10] + Qr.u[10] / 2 * (MPDT[i][j + 1].bz - MPDT[i][j].bz) + Qz.u[10] / 4 * (MPDT[i + 1][j].bz - 2 * MPDT[i][j].bz + MPDT[i - 1][j].bz);
					U[i][j].u[11] = U_bar2[i][j].u[11] + Qr.u[11] / 2 * (MPDT[i][j + 1].pe - MPDT[i][j].pe) + Qz.u[11] / 4 * (MPDT[i + 1][j].pe - 2 * MPDT[i][j].pe + MPDT[i - 1][j].pe);
					U[i][j].u[12] = U_bar2[i][j].u[12] + Qr.u[12] / 2 * (MPDT[i][j + 1].pi - MPDT[i][j].pi) + Qz.u[12] / 4 * (MPDT[i + 1][j].pi - 2 * MPDT[i][j].pi + MPDT[i - 1][j].pi);

				}
			}
			else //if (btype[i][j] == 0)
			{
				for (k = 0; k < 13; k++)
				{
					U[i][j].u[k] = U_bar2[i][j].u[k];
				}
			}
		}
	}

	//更新边界信息



	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			//if (U[i][j].u[0] > 0)
			//{
			//	MPDT[i][j].ne = U[i][j].u[0];
			//}
			//else
			//{
			//	MPDT[i][j].ne = 0;
			//}

			//if (U[i][j].u2 > 0)
			//{
			//	MPDT[i][j].ni = U[i][j].u2;
			//}
			//else
			//{
			//	MPDT[i][j].ni = 0;
			//}

			MPDT[i][j].ne = U[i][j].u[0];
			MPDT[i][j].ni = U[i][j].u[1];
			MPDT[i][j].ver = U[i][j].u[2];
			MPDT[i][j].vetheta = U[i][j].u[3];
			MPDT[i][j].vez = U[i][j].u[4];


			MPDT[i][j].vir = U[i][j].u[5];
			MPDT[i][j].vitheta = U[i][j].u[6];
			MPDT[i][j].viz = U[i][j].u[7];


			MPDT[i][j].br = U[i][j].u[8];
			MPDT[i][j].btheta = U[i][j].u[9];
			MPDT[i][j].bz = U[i][j].u[10];

			MPDT[i][j].pe = U[i][j].u[11];
			MPDT[i][j].pi = U[i][j].u[12];

			//MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) + 0.5 * MPDT[i][j].ne * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
			//MPDT[i][j].ei = MPDT[i][j].pi / (gamma - 1) + 0.5 * MPDT[i][j].ni * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);

		}
	}

	//smooth_ne(i, j);

	return;
}


//__global__ void cal_u_Kernel(int* c, const int* a, const int* b)
//{
//	int j = blockIdx.x * blockDim.x + threadIdx.x;
//	int i = threadIdx.x;
//	//c[i] = a[i] + b[i];
//	c[i] = j;
//}

struct _U cal_u(int i, int j)
{
	struct _U uij;
	uij.r = j;
	uij.z = i;
	uij.u[0] = MPDT[i][j].ne;
	uij.u[1] = MPDT[i][j].ni;
	uij.u[2] = MPDT[i][j].ver;
	uij.u[3] = MPDT[i][j].vetheta;
	uij.u[4] = MPDT[i][j].vez;
	uij.u[5] = MPDT[i][j].vir;
	uij.u[6] = MPDT[i][j].vitheta;
	uij.u[7] = MPDT[i][j].viz;
	uij.u[8] = MPDT[i][j].br;
	uij.u[9] = MPDT[i][j].btheta;
	uij.u[10] = MPDT[i][j].bz;
	uij.u[11] = MPDT[i][j].pe;
	uij.u[12] = MPDT[i][j].pi;
	return uij;
}

//__global__ void cal_fr_1_Kernel(const float * c, const int* a, const int* b)
//{
//	int j = blockIdx.x * blockDim.x + threadIdx.x;
//	int i = threadIdx.x;
//	//c[i] = a[i] + b[i];
//	c[i] = j;
//}

struct _F cal_fr(struct _U uij)
{
	struct _F fij;

	double ne = uij.u[0];
	double ni = uij.u[1];
	double ver = 0;
	double vetheta = 0;
	double vez = 0;
	double vir = 0;
	double vitheta = 0;
	double viz = 0;

	ver = uij.u[2];
	vetheta = uij.u[3];
	vez = uij.u[4];
	vir = uij.u[5];
	vitheta = uij.u[6];
	viz = uij.u[7];

	double br = uij.u[8];
	double btheta = uij.u[9];
	double bz = uij.u[10];

	double pe = uij.u[11];
	double pi = uij.u[12];

	//double ee = pe / (gamma - 1) + 0.5 * ne * ME * (ver * ver + vetheta * vetheta + vez * vez);
	//double ei = pi / (gamma - 1) + 0.5 * ni * MI * (vir * vir + vitheta * vitheta + viz * viz);

	fij.r = uij.r;
	fij.z = uij.z;
	fij.f[0] = ne * ver;
	fij.f[1] = ni * vir;
	if (ne == 0)
	{
		fij.f[2] = 0;
		fij.f[3] = 0;
		fij.f[4] = 0;
	}
	else
	{
		fij.f[2] = (pe + (btheta * btheta + bz * bz - br * br) / (2 * MU_0)) / (ne * ME);
		fij.f[3] = -(btheta * br) / (MU_0 * ne * ME);
		fij.f[4] = -(bz * br) / (MU_0 * ne * ME);
	}

	if (ni == 0)
	{
		fij.f[5] = 0;
		fij.f[6] = 0;
		fij.f[7] = 0;
	}
	else
	{
		fij.f[5] = (pi + (btheta * btheta + bz * bz - br * br) / (2 * MU_0)) / (ni * MI);
		fij.f[6] = -(btheta * br) / (MU_0 * ni * MI);
		fij.f[7] = -(bz * br) / (MU_0 * ni * MI);
	}

	fij.f[8] = 0;
	fij.f[9] = 0.5 * (vir * btheta - ver * btheta - br * vitheta + br * vetheta);
	fij.f[10] = 0.5 * (vir * bz - ver * bz - br * viz + br * vez);
	fij.f[11] = pe * ver + (gamma - 1) * (-ver * (btheta * btheta + bz * bz - br * br) / (2 * MU_0) + vetheta * btheta * br / MU_0 + vez * bz * br / MU_0);
	fij.f[12] = pi * vir + (gamma - 1) * (-vir * (btheta * btheta + bz * bz - br * br) / (2 * MU_0) + vitheta * btheta * br / MU_0 + viz * bz * br / MU_0);
	return fij;
}

struct _F cal_fz(struct _U uij)
{
	struct _F fij;

	double ne = uij.u[0];
	double ni = uij.u[1];
	double ver = 0;
	double vetheta = 0;
	double vez = 0;
	double vir = 0;
	double vitheta = 0;
	double viz = 0;

	ver = uij.u[2];
	vetheta = uij.u[3];
	vez = uij.u[4];
	vir = uij.u[5];
	vitheta = uij.u[6];
	viz = uij.u[7];

	double br = uij.u[8];
	double btheta = uij.u[9];
	double bz = uij.u[10];

	double pe = uij.u[11];
	double pi = uij.u[12];

	//double ee = pe / (gamma - 1) + 0.5 * ne * ME * (ver * ver + vetheta * vetheta + vez * vez);
	//double ei = pi / (gamma - 1) + 0.5 * ni * MI * (vir * vir + vitheta * vitheta + viz * viz);


	fij.r = uij.r;
	fij.z = uij.z;

	fij.f[0] = ne * vez;
	fij.f[1] = ni * viz;
	if (ne == 0)
	{
		fij.f[2] = 0;
		fij.f[3] = 0;
		fij.f[4] = 0;
	}
	else
	{
		fij.f[2] = -(bz * br) / (MU_0 * ne * ME);
		fij.f[3] = -(btheta * bz) / (MU_0 * ne * ME);
		fij.f[4] = (pe + (btheta * btheta + br * br - bz * bz) / (2 * MU_0)) / (ne * ME);
	}

	if (ni == 0)
	{
		fij.f[5] = 0;
		fij.f[6] = 0;
		fij.f[7] = 0;
	}
	else
	{
		fij.f[5] = -(bz * br) / (MU_0 * ni * MI);
		fij.f[6] = -(btheta * bz) / (MU_0 * ni * MI);
		fij.f[7] = (pi + (btheta * btheta + br * br - bz * bz) / (2 * MU_0)) / (ni * MI);
	}

	fij.f[8] = 0.5 * (viz * br - vez * br + bz * ver - bz * vir);
	fij.f[9] = 0.5 * (viz * btheta - vez * btheta - bz * vitheta + bz * vetheta);
	fij.f[10] = 0;
	fij.f[11] = pe * vez + (gamma - 1) * (-vez * (btheta * btheta + br * br - bz * bz) / (2 * MU_0) + vetheta * btheta * bz / MU_0 + ver * bz * br / MU_0);
	fij.f[12] = pi * viz + (gamma - 1) * (-viz * (btheta * btheta + br * br - bz * bz) / (2 * MU_0) + vitheta * btheta * bz / MU_0 + vir * bz * br / MU_0);
	return fij;
}


struct _F cal_s(struct _U uij)
{
	struct _F fij;

	double ne = uij.u[0];
	double ni = uij.u[1];
	double ver = 0;
	double vetheta = 0;
	double vez = 0;
	double vir = 0;
	double vitheta = 0;
	double viz = 0;

	ver = uij.u[2];
	vetheta = uij.u[3];
	vez = uij.u[4];
	vir = uij.u[5];
	vitheta = uij.u[6];
	viz = uij.u[7];

	double br = uij.u[8];
	double btheta = uij.u[9];
	double bz = uij.u[10];

	double pe = uij.u[11];
	double pi = uij.u[12];

	double ee = pe / (gamma - 1) + 0.5 * ne * ME * (ver * ver + vetheta * vetheta + vez * vez);
	double ei = pi / (gamma - 1) + 0.5 * ni * MI * (vir * vir + vitheta * vitheta + viz * viz);

	fij.r = uij.r;
	fij.z = uij.z;
	double r1 = 0;
	if (uij.r == 0)
	{
		r1 = -2 / dr;
	}
	else 
	{
		r1 = -1 / (uij.r * dr);
	}
	fij.f[0] = r1 * ne * ver;
	fij.f[1] = r1 * ni * vir;

	if (ne == 0)
	{
		fij.f[2] = 0;
		fij.f[3] = 0;
		fij.f[4] = 0;
	}
	else
	{
		fij.f[2] = r1 * (pe + (btheta * btheta + bz * bz - br * br) / (2 * MU_0)) / (ne * ME);
		fij.f[3] = -r1 * (btheta * br) / (MU_0 * ne * ME);
		fij.f[4] = -r1 * (bz * br) / (MU_0 * ne * ME);
	}

	if (ni == 0)
	{
		fij.f[5] = 0;
		fij.f[6] = 0;
		fij.f[7] = 0;
	}
	else
	{
		fij.f[5] = r1 * (pi + (btheta * btheta + bz * bz - br * br) / (2 * MU_0)) / (ni * MI);
		fij.f[6] = -r1 * (btheta * br) / (MU_0 * ni * MI);
		fij.f[7] = -r1 * (bz * br) / (MU_0 * ni * MI);
	}

	fij.f[8] = 0;
	fij.f[9] = r1 * 0.5 * (vir * btheta - ver * btheta - br * vitheta + br * vetheta);
	fij.f[10] = r1 * 0.5 * (vir * bz - ver * bz - br * viz + br * vez);
	fij.f[11] = r1 * (pe * ver + (gamma - 1) * (-ver * (btheta * btheta + bz * bz - br * br) / (2 * MU_0) + vetheta * btheta * br / MU_0 + vez * bz * br / MU_0));
	fij.f[12] = r1 * (pi * vir + (gamma - 1) * (-vir * (btheta * btheta + bz * bz - br * br) / (2 * MU_0) + vitheta * btheta * br / MU_0 + viz * bz * br / MU_0));
	return fij;
}

struct _U arti_vis(struct  node np, struct  node n, struct  node nn)
{
	struct _U Q;
	double eta = 0.5;

	for (int k = 0; k < 13; k++)
	{
		Q.u[k] = 0;
	}

	if (abs(np.ne - n.ne) + abs(n.ne - nn.ne) != 0)
	{
		Q.u[0] = eta * abs(abs(np.ne - n.ne) - abs(n.ne - nn.ne)) / abs(abs(np.ne - n.ne) + abs(n.ne - nn.ne));
	}

	if (abs(np.ni - n.ni) + abs(n.ni - nn.ni) != 0)
	{
		Q.u[1] = eta * abs(abs(np.ni - n.ni) - abs(n.ni - nn.ni)) / abs(abs(np.ni - n.ni) + abs(n.ni - nn.ni));
	}

	Q.u[2] = Q.u[0];
	Q.u[3] = Q.u[0];
	Q.u[4] = Q.u[0];
	Q.u[11] = Q.u[0];

	Q.u[5] = Q.u[1];
	Q.u[6] = Q.u[1];
	Q.u[7] = Q.u[1];
	Q.u[12] = Q.u[1];

	return Q;
}

//struct _U arti_vis(struct  node np, struct  node n, struct  node nn)
//{
//	struct _U Q;
//	double eta = 0.5;
//
//	for (int k = 0; k < 13; k++)
//	{
//		Q.u[k] = 0;
//	}
//
//	if (abs(np.ne - n.ne) + abs(n.ne - nn.ne) != 0)
//	{
//		Q.u1 = eta * abs(abs(np.ne - n.ne) - abs(n.ne - nn.ne)) / abs(abs(np.ne - n.ne) + abs(n.ne - nn.ne));
//	}
//
//	if (abs(np.ni - n.ni) + abs(n.ni - nn.ni) != 0)
//	{
//		Q.u2 = eta * abs(abs(np.ni - n.ni) - abs(n.ni - nn.ni)) / abs(abs(np.ni - n.ni) + abs(n.ni - nn.ni));
//	}
//
//	if (abs(np.ver - n.ver) + abs(n.ver - nn.ver) != 0)
//	{
//		Q.u3 = eta * abs(abs(np.ver - n.ver) - abs(n.ver - nn.ver)) / abs(abs(np.ver - n.ver) + abs(n.ver - nn.ver));
//	}
//
//	if (abs(np.vetheta - n.vetheta) + abs(n.vetheta - nn.vetheta) != 0)
//	{
//		Q.u4 = eta * abs(abs(np.vetheta - n.vetheta) - abs(n.vetheta - nn.vetheta)) / abs(abs(np.vetheta - n.vetheta) + abs(n.vetheta - nn.vetheta));
//	}
//
//	if (abs(np.vez - n.vez) + abs(n.vez - nn.vez) != 0)
//	{
//		Q.u5 = eta * abs(abs(np.vez - n.vez) - abs(n.vez - nn.vez)) / abs(abs(np.vez - n.vez) + abs(n.vez - nn.vez));
//	}
//
//	if (abs(np.vir - n.vir) + abs(n.vir - nn.vir) != 0)
//	{
//		Q.u6 = eta * abs(abs(np.vir - n.vir) - abs(n.vir - nn.vir)) / abs(abs(np.vir - n.vir) + abs(n.vir - nn.vir));
//	}
//
//	if (abs(np.vitheta - n.vitheta) + abs(n.vitheta - nn.vitheta) != 0)
//	{
//		Q.u7 = eta * abs(abs(np.vitheta - n.vitheta) - abs(n.vitheta - nn.vitheta)) / abs(abs(np.vitheta - n.vitheta) + abs(n.vitheta - nn.vitheta));
//	}
//
//	if (abs(np.viz - n.viz) + abs(n.viz - nn.viz) != 0)
//	{
//		Q.u8 = eta * abs(abs(np.viz - n.viz) - abs(n.viz - nn.viz)) / abs(abs(np.viz - n.viz) + abs(n.viz - nn.viz));
//	}
//
//	return Q;
//}

//计算多余电荷
//void Q_fluid()
//{
//	struct _U Qr;
//	struct _U Qz;
//	register int i, j, k;
//	for (i = 0; i < nz; i++)
//	{
//		for (j = 0; j < nr; j++)
//		{
//
//			Uq[i][j] = cal_qu(i, j);
//			Fqr[i][j] = cal_fqr(Uq[i][j]);
//			Fqz[i][j] = cal_fqz(Uq[i][j]);
//			sq[i][j] = cal_qs(Uq[i][j]);
//
//			Uq_bar[i][j].z = i;
//			Uq_bar[i][j].r = j;
//			Uq_bar2[i][j].z = i;
//			Uq_bar2[i][j].r = j;
//
//		}
//	}
//
//	for (i = 0; i < nz; i++)
//	{
//		for (j = 0; j < nr; j++)
//		{
//			if (btype[i][j] == 1)
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar[i][j].u[k] = Uq[i][j].u[k] - dt / dr * (Fqr[i][j + 1].f[k] - Fqr[i][j].f[k]) - dt / dz * (Fqz[i + 1][j].f[k] - Fqz[i][j].f[k]) + dt * sq[i][j].f[k];
//				}
//			}
//			else if (btype[i][j] == LEFT)//左边的边界，复制右边的参数
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar[i][j].u[k] = Uq[i][j].u[k] - dt / dr * (Fqr[i][j + 1].f[k] - Fqr[i][j].f[k]) - dt / dz * (Fqz[i + 1][j].f[k] - Fqz[i][j].f[k]) + dt * sq[i][j].f[k];
//				}
//			}
//			else if (btype[i][j] == RIGHT)
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar[i][j].u[k] = Uq[i][j].u[k] - dt / dr * (Fqr[i][j + 1].f[k] - Fqr[i][j].f[k]) - dt / dz * (Fqz[i][j].f[k] - Fqz[i - 1][j].f[k]) + dt * sq[i][j].f[k];
//				}
//			}
//			else if (btype[i][j] == UP)
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar[i][j].u[k] = Uq[i][j].u[k] - dt / dr * (Fqr[i][j].f[k] - Fqr[i][j - 1].f[k]) - dt / dz * (Fqz[i + 1][j].f[k] - Fqz[i][j].f[k]) + dt * sq[i][j].f[k];
//				}
//			}
//			else if (btype[i][j] == DOWN)
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar[i][j].u[k] = Uq[i][j].u[k] - dt / dr * (Fqr[i][j + 1].f[k] - Fqr[i][j].f[k]) - dt / dz * (Fqz[i + 1][j].f[k] - Fqz[i][j].f[k]) + dt * sq[i][j].f[k];
//				}
//			}
//			else if (btype[i][j] == (LEFT + UP))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar[i][j].u[k] = Uq[i][j].u[k] - dt / dr * (Fqr[i][j].f[k] - Fqr[i][j - 1].f[k]) - dt / dz * (Fqz[i + 1][j].f[k] - Fqz[i][j].f[k]) + dt * sq[i][j].f[k];
//				}
//			}
//			else if (btype[i][j] == (LEFT + DOWN))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar[i][j].u[k] = Uq[i][j].u[k] - dt / dr * (Fqr[i][j + 1].f[k] - Fqr[i][j].f[k]) - dt / dz * (Fqz[i + 1][j].f[k] - Fqz[i][j].f[k]) + dt * sq[i][j].f[k];
//				}
//			}
//			else if (btype[i][j] == (RIGHT + DOWN))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar[i][j].u[k] = Uq[i][j].u[k] - dt / dr * (Fqr[i][j + 1].f[k] - Fqr[i][j].f[k]) - dt / dz * (Fqz[i][j].f[k] - Fqz[i - 1][j].f[k]) + dt * sq[i][j].f[k];
//				}
//			}
//			else if (btype[i][j] == (RIGHT + UP))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar[i][j].u[k] = Uq[i][j].u[k] - dt / dr * (Fqr[i][j].f[k] - Fqr[i][j - 1].f[k]) - dt / dz * (Fqz[i][j].f[k] - Fqz[i - 1][j].f[k]) + dt * sq[i][j].f[k];
//				}
//			}
//			else if (btype[i][j] == 0)
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar[i][j].u[k] = 0;
//				}
//			}
//
//			Fqr_bar[i][j] = cal_fqr(Uq_bar[i][j]);
//			Fqz_bar[i][j] = cal_fqz(Uq_bar[i][j]);
//			sq_bar[i][j] = cal_qs(Uq_bar[i][j]);
//		}
//	}
//
//	//矫正步
//	for (i = 0; i < nz; i++)
//	{
//		for (j = 0; j < nr; j++)
//		{
//			if (btype[i][j] == 1)
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar2[i][j].u[k] = 0.5 * (Uq[i][j].u[k] + Uq_bar[i][j].u[k] - dt / dr * (Fqr_bar[i][j].f[k] - Fqr_bar[i][j - 1].f[k]) - dt / dz * (Fqz_bar[i][j].f[k] - Fqz_bar[i - 1][j].f[k]) + dt * sq_bar[i][j].f[k]);
//				}
//			}
//			else if (btype[i][j] == LEFT)//左侧边界
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar2[i][j].u[k] = 0.5 * (Uq[i][j].u[k] + Uq_bar[i][j].u[k] - dt / dr * (Fqr_bar[i][j].f[k] - Fqr_bar[i][j - 1].f[k]) - dt / dz * (Fqz_bar[i + 1][j].f[k] - Fqz_bar[i][j].f[k]) + dt * sq_bar[i][j].f[k]);
//				}
//			}
//			else if (btype[i][j] == RIGHT)
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar2[i][j].u[k] = 0.5 * (Uq[i][j].u[k] + Uq_bar[i][j].u[k] - dt / dr * (Fqr_bar[i][j].f[k] - Fqr_bar[i][j - 1].f[k]) - dt / dz * (Fqz_bar[i][j].f[k] - Fqz_bar[i - 1][j].f[k]) + dt * sq_bar[i][j].f[k]);
//				}
//			}
//			else if (btype[i][j] == UP)
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar2[i][j].u[k] = 0.5 * (Uq[i][j].u[k] + Uq_bar[i][j].u[k] - dt / dr * (Fqr_bar[i][j].f[k] - Fqr_bar[i][j - 1].f[k]) - dt / dz * (Fqz_bar[i][j].f[k] - Fqz_bar[i - 1][j].f[k]) + dt * sq_bar[i][j].f[k]);
//				}
//			}
//			else if (btype[i][j] == DOWN)
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar2[i][j].u[k] = 0.5 * (Uq[i][j].u[k] + Uq_bar[i][j].u[k] - dt / dr * (Fqr_bar[i][j + 1].f[k] - Fqr_bar[i][j].f[k]) - dt / dz * (Fqz_bar[i][j].f[k] - Fqz_bar[i - 1][j].f[k]) + dt * sq_bar[i][j].f[k]);
//				}
//			}
//			else if (btype[i][j] == (LEFT + UP))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar2[i][j].u[k] = 0.5 * (Uq[i][j].u[k] + Uq_bar[i][j].u[k] - dt / dr * (Fqr_bar[i][j].f[k] - Fqr_bar[i][j - 1].f[k]) - dt / dz * (Fqz_bar[i + 1][j].f[k] - Fqz_bar[i][j].f[k]) + dt * sq_bar[i][j].f[k]);
//				}
//			}
//			else if (btype[i][j] == (LEFT + DOWN))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar2[i][j].u[k] = 0.5 * (Uq[i][j].u[k] + Uq_bar[i][j].u[k] - dt / dr * (Fqr_bar[i][j + 1].f[k] - Fqr_bar[i][j].f[k]) - dt / dz * (Fqz_bar[i + 1][j].f[k] - Fqz_bar[i][j].f[k]) + dt * sq_bar[i][j].f[k]);
//				}
//			}
//			else if (btype[i][j] == (RIGHT + DOWN))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar2[i][j].u[k] = 0.5 * (Uq[i][j].u[k] + Uq_bar[i][j].u[k] - dt / dr * (Fqr_bar[i][j + 1].f[k] - Fqr_bar[i][j].f[k]) - dt / dz * (Fqz_bar[i][j].f[k] - Fqz_bar[i - 1][j].f[k]) + dt * sq_bar[i][j].f[k]);
//				}
//			}
//			else if (btype[i][j] == (RIGHT + UP))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar2[i][j].u[k] = 0.5 * (Uq[i][j].u[k] + Uq_bar[i][j].u[k] - dt / dr * (Fqr_bar[i][j].f[k] - Fqr_bar[i][j - 1].f[k]) - dt / dz * (Fqz_bar[i][j].f[k] - Fqz_bar[i - 1][j].f[k]) + dt * sq_bar[i][j].f[k]);
//				}
//			}
//			else if (btype[i][j] == 0)
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq_bar2[i][j].u[k] = 0;
//				}
//			}
//
//		}
//	}
//
//	//计算				
//
//	for (i = 0; i < nz; i++)
//	{
//		for (j = 0; j < nr; j++)
//		{
//			if (btype[i][j] == 1)
//			{
//				Qz = arti_q_vis(MPDT[i + 1][j], MPDT[i][j], MPDT[i - 1][j]);
//				Qr = arti_q_vis(MPDT[i][j + 1], MPDT[i][j], MPDT[i][j - 1]);
//
//				Uq[i][j].u[0] = Uq_bar2[i][j].u[0] + Qr.u[0] / 4 * (MPDT[i][j + 1].neq - 2 * MPDT[i][j].neq + MPDT[i][j - 1].neq) + Qz.u[0] / 4 * (MPDT[i + 1][j].neq - 2 * MPDT[i][j].neq + MPDT[i - 1][j].neq);
//				Uq[i][j].u[1] = Uq_bar2[i][j].u[1] + Qr.u[1] / 4 * (MPDT[i][j + 1].peq - 2 * MPDT[i][j].peq + MPDT[i][j - 1].peq) + Qz.u[1] / 4 * (MPDT[i + 1][j].peq - 2 * MPDT[i][j].peq + MPDT[i - 1][j].peq);
//				Uq[i][j].u[2] = Uq_bar2[i][j].u[2] + Qr.u[2] / 4 * (MPDT[i][j + 1].vnqr - 2 * MPDT[i][j].vnqr + MPDT[i][j - 1].vnqr) + Qz.u[2] / 4 * (MPDT[i + 1][j].vnqr - 2 * MPDT[i][j].vnqr + MPDT[i - 1][j].vnqr);
//				Uq[i][j].u[3] = Uq_bar2[i][j].u[3] + Qr.u[3] / 4 * (MPDT[i][j + 1].vnqtheta - 2 * MPDT[i][j].vnqtheta + MPDT[i][j - 1].vnqtheta) + Qz.u[3] / 4 * (MPDT[i + 1][j].vnqtheta - 2 * MPDT[i][j].vnqtheta + MPDT[i - 1][j].vnqtheta);
//				Uq[i][j].u[4] = Uq_bar2[i][j].u[4] + Qr.u[4] / 4 * (MPDT[i][j + 1].vnqz - 2 * MPDT[i][j].vnqz + MPDT[i][j - 1].vnqz) + Qz.u[4] / 4 * (MPDT[i + 1][j].vnqz - 2 * MPDT[i][j].vnqz + MPDT[i - 1][j].vnqz);
//				Uq[i][j].u[5] = Uq_bar2[i][j].u[5] + Qr.u[5] / 4 * (MPDT[i][j + 1].vpqr - 2 * MPDT[i][j].vpqr + MPDT[i][j - 1].vpqr) + Qz.u[5] / 4 * (MPDT[i + 1][j].vpqr - 2 * MPDT[i][j].vpqr + MPDT[i - 1][j].vpqr);
//				Uq[i][j].u[6] = Uq_bar2[i][j].u[6] + Qr.u[6] / 4 * (MPDT[i][j + 1].vpqtheta - 2 * MPDT[i][j].vpqtheta + MPDT[i][j - 1].vpqtheta) + Qz.u[6] / 4 * (MPDT[i + 1][j].vpqtheta - 2 * MPDT[i][j].vpqtheta + MPDT[i - 1][j].vpqtheta);
//				Uq[i][j].u[7] = Uq_bar2[i][j].u[7] + Qr.u[7] / 4 * (MPDT[i][j + 1].vpqz - 2 * MPDT[i][j].vpqz + MPDT[i][j - 1].vpqz) + Qz.u[7] / 4 * (MPDT[i + 1][j].vpqz - 2 * MPDT[i][j].vpqz + MPDT[i - 1][j].vpqz);
//
//
//#ifdef FLUID_DEBUG
//				if (i == row && j == col)
//				{
//					//cout << "final U[i][j].u1 = " << U[i][j].u1 << endl;
//					printf("final U[i][j].u1 = %.5e\n", U[i][j].u[0]);
//					printf("U_bar2[i][j + 1].u1 = %.5e\n", U_bar2[i][j + 1].u[0]);
//					printf("U_bar2[i][j - 1].u1 = %.5e\n", U_bar2[i][j - 1].u[0]);
//					printf("U_bar2[i + 1][j].u1 = %.5e\n", U_bar2[i + 1][j].u[0]);
//					printf("U_bar2[i - 1][j].u1 = %.5e\n", U_bar2[i - 1][j].u[0]);
//
//					printf("final U[i][j].u3 = %.5e\n", U[i][j].u[2]);
//					printf("U_bar2[i][j + 1].u3 = %.5e\n", U_bar2[i][j + 1].u[2]);
//					printf("U_bar2[i][j - 1].u3 = %.5e\n", U_bar2[i][j - 1].u[2]);
//					printf("U_bar2[i + 1][j].u3 = %.5e\n", U_bar2[i + 1][j].u[2]);
//					printf("U_bar2[i - 1][j].u3 = %.5e\n", U_bar2[i - 1][j].u[2]);
//
//				}
//#endif// FLUID_DEBUG	
//
//			}
//			else if (btype[i][j] == LEFT)//左边的边界，复制右边的参数
//			{
//				Qr = arti_q_vis(MPDT[i][j + 1], MPDT[i][j], MPDT[i][j - 1]);
//
//
//				Uq[i][j].u[0] = Uq_bar2[i][j].u[0] + Qr.u[0] / 4 * (MPDT[i][j + 1].neq - 2 * MPDT[i][j].neq + MPDT[i][j - 1].neq);
//				Uq[i][j].u[1] = Uq_bar2[i][j].u[1] + Qr.u[1] / 4 * (MPDT[i][j + 1].peq - 2 * MPDT[i][j].peq + MPDT[i][j - 1].peq);
//				Uq[i][j].u[2] = Uq_bar2[i][j].u[2] + Qr.u[2] / 4 * (MPDT[i][j + 1].vnqr - 2 * MPDT[i][j].vnqr + MPDT[i][j - 1].vnqr);
//				Uq[i][j].u[3] = Uq_bar2[i][j].u[3] + Qr.u[3] / 4 * (MPDT[i][j + 1].vnqtheta - 2 * MPDT[i][j].vnqtheta + MPDT[i][j - 1].vnqtheta);
//				Uq[i][j].u[4] = Uq_bar2[i][j].u[4] + Qr.u[4] / 4 * (MPDT[i][j + 1].vnqz - 2 * MPDT[i][j].vnqz + MPDT[i][j - 1].vnqz);
//				Uq[i][j].u[5] = Uq_bar2[i][j].u[5] + Qr.u[5] / 4 * (MPDT[i][j + 1].vpqr - 2 * MPDT[i][j].vpqr + MPDT[i][j - 1].vpqr);
//				Uq[i][j].u[6] = Uq_bar2[i][j].u[6] + Qr.u[6] / 4 * (MPDT[i][j + 1].vpqtheta - 2 * MPDT[i][j].vpqtheta + MPDT[i][j - 1].vpqtheta);
//				Uq[i][j].u[7] = Uq_bar2[i][j].u[7] + Qr.u[7] / 4 * (MPDT[i][j + 1].vpqz - 2 * MPDT[i][j].vpqz + MPDT[i][j - 1].vpqz);
//
//			}
//			else if (btype[i][j] == RIGHT)
//			{
//
//				Qr = arti_q_vis(MPDT[i][j + 1], MPDT[i][j], MPDT[i][j - 1]);
//
//
//				Uq[i][j].u[0] = Uq_bar2[i][j].u[0] + Qr.u[0] / 4 * (MPDT[i][j + 1].neq - 2 * MPDT[i][j].neq + MPDT[i][j - 1].neq);
//				Uq[i][j].u[1] = Uq_bar2[i][j].u[1] + Qr.u[1] / 4 * (MPDT[i][j + 1].peq - 2 * MPDT[i][j].peq + MPDT[i][j - 1].peq);
//				Uq[i][j].u[2] = Uq_bar2[i][j].u[2] + Qr.u[2] / 4 * (MPDT[i][j + 1].vnqr - 2 * MPDT[i][j].vnqr + MPDT[i][j - 1].vnqr);
//				Uq[i][j].u[3] = Uq_bar2[i][j].u[3] + Qr.u[3] / 4 * (MPDT[i][j + 1].vnqtheta - 2 * MPDT[i][j].vnqtheta + MPDT[i][j - 1].vnqtheta);
//				Uq[i][j].u[4] = Uq_bar2[i][j].u[4] + Qr.u[4] / 4 * (MPDT[i][j + 1].vnqz - 2 * MPDT[i][j].vnqz + MPDT[i][j - 1].vnqz);
//				Uq[i][j].u[5] = Uq_bar2[i][j].u[5] + Qr.u[5] / 4 * (MPDT[i][j + 1].vpqr - 2 * MPDT[i][j].vpqr + MPDT[i][j - 1].vpqr);
//				Uq[i][j].u[6] = Uq_bar2[i][j].u[6] + Qr.u[6] / 4 * (MPDT[i][j + 1].vpqtheta - 2 * MPDT[i][j].vpqtheta + MPDT[i][j - 1].vpqtheta);
//				Uq[i][j].u[7] = Uq_bar2[i][j].u[7] + Qr.u[7] / 4 * (MPDT[i][j + 1].vpqz - 2 * MPDT[i][j].vpqz + MPDT[i][j - 1].vpqz);
//			}
//			else if (btype[i][j] == UP)
//			{
//				Qz = arti_q_vis(MPDT[i + 1][j], MPDT[i][j], MPDT[i - 1][j]);
//
//				Uq[i][j].u[0] = Uq_bar2[i][j].u[0] + Qz.u[0] / 4 * (MPDT[i][j + 1].neq - 2 * MPDT[i][j].neq + MPDT[i][j - 1].neq);
//				Uq[i][j].u[1] = Uq_bar2[i][j].u[1] + Qz.u[1] / 4 * (MPDT[i][j + 1].peq - 2 * MPDT[i][j].peq + MPDT[i][j - 1].peq);
//				Uq[i][j].u[2] = Uq_bar2[i][j].u[2] + Qz.u[2] / 4 * (MPDT[i][j + 1].vnqr - 2 * MPDT[i][j].vnqr + MPDT[i][j - 1].vnqr);
//				Uq[i][j].u[3] = Uq_bar2[i][j].u[3] + Qz.u[3] / 4 * (MPDT[i][j + 1].vnqtheta - 2 * MPDT[i][j].vnqtheta + MPDT[i][j - 1].vnqtheta);
//				Uq[i][j].u[4] = Uq_bar2[i][j].u[4] + Qz.u[4] / 4 * (MPDT[i][j + 1].vnqz - 2 * MPDT[i][j].vnqz + MPDT[i][j - 1].vnqz);
//				Uq[i][j].u[5] = Uq_bar2[i][j].u[5] + Qz.u[5] / 4 * (MPDT[i][j + 1].vpqr - 2 * MPDT[i][j].vpqr + MPDT[i][j - 1].vpqr);
//				Uq[i][j].u[6] = Uq_bar2[i][j].u[6] + Qz.u[6] / 4 * (MPDT[i][j + 1].vpqtheta - 2 * MPDT[i][j].vpqtheta + MPDT[i][j - 1].vpqtheta);
//				Uq[i][j].u[7] = Uq_bar2[i][j].u[7] + Qz.u[7] / 4 * (MPDT[i][j + 1].vpqz - 2 * MPDT[i][j].vpqz + MPDT[i][j - 1].vpqz);
//
//
//
//			}
//			else if (btype[i][j] == DOWN)
//			{
//				Qz = arti_q_vis(MPDT[i + 1][j], MPDT[i][j], MPDT[i - 1][j]);
//
//				Uq[i][j].u[0] = Uq_bar2[i][j].u[0] + Qz.u[0] / 4 * (MPDT[i][j + 1].neq - 2 * MPDT[i][j].neq + MPDT[i][j - 1].neq);
//				Uq[i][j].u[1] = Uq_bar2[i][j].u[1] + Qz.u[1] / 4 * (MPDT[i][j + 1].peq - 2 * MPDT[i][j].peq + MPDT[i][j - 1].peq);
//				Uq[i][j].u[2] = Uq_bar2[i][j].u[2] + Qz.u[2] / 4 * (MPDT[i][j + 1].vnqr - 2 * MPDT[i][j].vnqr + MPDT[i][j - 1].vnqr);
//				Uq[i][j].u[3] = Uq_bar2[i][j].u[3] + Qz.u[3] / 4 * (MPDT[i][j + 1].vnqtheta - 2 * MPDT[i][j].vnqtheta + MPDT[i][j - 1].vnqtheta);
//				Uq[i][j].u[4] = Uq_bar2[i][j].u[4] + Qz.u[4] / 4 * (MPDT[i][j + 1].vnqz - 2 * MPDT[i][j].vnqz + MPDT[i][j - 1].vnqz);
//				Uq[i][j].u[5] = Uq_bar2[i][j].u[5] + Qz.u[5] / 4 * (MPDT[i][j + 1].vpqr - 2 * MPDT[i][j].vpqr + MPDT[i][j - 1].vpqr);
//				Uq[i][j].u[6] = Uq_bar2[i][j].u[6] + Qz.u[6] / 4 * (MPDT[i][j + 1].vpqtheta - 2 * MPDT[i][j].vpqtheta + MPDT[i][j - 1].vpqtheta);
//				Uq[i][j].u[7] = Uq_bar2[i][j].u[7] + Qz.u[7] / 4 * (MPDT[i][j + 1].vpqz - 2 * MPDT[i][j].vpqz + MPDT[i][j - 1].vpqz);
//			}
//			else if (btype[i][j] == (LEFT + UP))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq[i][j].u[k] = Uq_bar2[i][j].u[k];
//				}
//			}
//			else if (btype[i][j] == (LEFT + DOWN))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq[i][j].u[k] = Uq_bar2[i][j].u[k];
//				}
//			}
//			else if (btype[i][j] == (RIGHT + DOWN))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq[i][j].u[k] = Uq_bar2[i][j].u[k];
//				}
//			}
//			else if (btype[i][j] == (RIGHT + UP))
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq[i][j].u[k] = Uq_bar2[i][j].u[k];
//				}
//			}
//			else //if (btype[i][j] == 0)
//			{
//				for (k = 0; k < 13; k++)
//				{
//					Uq[i][j].u[k] = Uq_bar2[i][j].u[k];
//				}
//			}
//		}
//	}
//
//	//更新边界信息
//
//
//
//	for (i = 0; i < nz; i++)
//	{
//		for (j = 0; j < nr; j++)
//		{
//			//if (U[i][j].u[0] > 0)
//			//{
//			//	MPDT[i][j].ne = U[i][j].u[0];
//			//}
//			//else
//			//{
//			//	MPDT[i][j].ne = 0;
//			//}
//
//			//if (U[i][j].u2 > 0)
//			//{
//			//	MPDT[i][j].ni = U[i][j].u2;
//			//}
//			//else
//			//{
//			//	MPDT[i][j].ni = 0;
//			//}
//
//			MPDT[i][j].neq = Uq[i][j].u[0];
//			MPDT[i][j].peq = Uq[i][j].u[1];
//			MPDT[i][j].vnqr = Uq[i][j].u[2];
//			MPDT[i][j].vnqtheta = Uq[i][j].u[3];
//			MPDT[i][j].vnqz = Uq[i][j].u[4];
//
//
//			MPDT[i][j].vpqr = Uq[i][j].u[5];
//			MPDT[i][j].vpqtheta = Uq[i][j].u[6];
//			MPDT[i][j].vpqz = Uq[i][j].u[7];
//
//
//		}
//	}
//
//
//	return;
//
//}

void Q_fluid()
{
	struct _U Qr;
	struct _U Qz;
	register int i, j, k;
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{

			Uq[i][j] = cal_qu(i, j);
			Fqr[i][j] = cal_fqr(Uq[i][j]);
			Fqz[i][j] = cal_fqz(Uq[i][j]);
			sq[i][j] = cal_qs(Uq[i][j]);

			Uq_bar[i][j].z = i;
			Uq_bar[i][j].r = j;
			Uq_bar2[i][j].z = i;
			Uq_bar2[i][j].r = j;

		}
	}

	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				for (k = 0; k < 13; k++)
				{
					Uq_bar[i][j].u[k] = Uq[i][j].u[k] - dt / dr * (Fqr[i][j + 1].f[k] - Fqr[i][j].f[k]) - dt / dz * (Fqz[i + 1][j].f[k] - Fqz[i][j].f[k]) + dt * sq[i][j].f[k];
				}
			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				for (k = 0; k < 13; k++)
				{
					Uq_bar[i][j].u[k] = Uq[i][j].u[k] - dt / dz * (Fqz[i + 1][j].f[k] - Fqz[i][j].f[k]) + dt * sq[i][j].f[k];
				}
			}
			else //if (btype[i][j] == 0)
			{
				for (k = 0; k < 13; k++)
				{
					Uq_bar[i][j].u[k] = Uq[i][j].u[k];
				}
			}

			Fqr_bar[i][j] = cal_fqr(Uq_bar[i][j]);
			Fqz_bar[i][j] = cal_fqz(Uq_bar[i][j]);
			sq_bar[i][j] = cal_qs(Uq_bar[i][j]);
		}
	}

	//矫正步
	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				for (k = 0; k < 13; k++)
				{
					Uq_bar2[i][j].u[k] = 0.5 * (Uq[i][j].u[k] + Uq_bar[i][j].u[k] - dt / dr * (Fqr_bar[i][j].f[k] - Fqr_bar[i][j - 1].f[k]) - dt / dz * (Fqz_bar[i][j].f[k] - Fqz_bar[i - 1][j].f[k]) + dt * sq_bar[i][j].f[k]);
				}
			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				for (k = 0; k < 13; k++)
				{
					Uq_bar2[i][j].u[k] = 0.5 * (Uq[i][j].u[k] + Uq_bar[i][j].u[k] -  dt / dz * (Fqz_bar[i][j].f[k] - Fqz_bar[i - 1][j].f[k]) + dt * sq_bar[i][j].f[k]);
				}
			}
			else //if (btype[i][j] == 0)
			{
				for (k = 0; k < 13; k++)
				{
					Uq_bar2[i][j].u[k] = Uq_bar[i][j].u[k];
				}
			}
		}
	}

	//计算				

	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (btype[i][j] == 1)
			{
				Qz = arti_q_vis(MPDT[i + 1][j], MPDT[i][j], MPDT[i - 1][j]);
				Qr = arti_q_vis(MPDT[i][j + 1], MPDT[i][j], MPDT[i][j - 1]);

				Uq[i][j].u[0] = Uq_bar2[i][j].u[0] + Qr.u[0] / 4 * (MPDT[i][j + 1].neq - 2 * MPDT[i][j].neq + MPDT[i][j - 1].neq) + Qz.u[0] / 4 * (MPDT[i + 1][j].neq - 2 * MPDT[i][j].neq + MPDT[i - 1][j].neq);
				Uq[i][j].u[1] = Uq_bar2[i][j].u[1] + Qr.u[1] / 4 * (MPDT[i][j + 1].peq - 2 * MPDT[i][j].peq + MPDT[i][j - 1].peq) + Qz.u[1] / 4 * (MPDT[i + 1][j].peq - 2 * MPDT[i][j].peq + MPDT[i - 1][j].peq);
				Uq[i][j].u[2] = Uq_bar2[i][j].u[2] + Qr.u[2] / 4 * (MPDT[i][j + 1].vnqr - 2 * MPDT[i][j].vnqr + MPDT[i][j - 1].vnqr) + Qz.u[2] / 4 * (MPDT[i + 1][j].vnqr - 2 * MPDT[i][j].vnqr + MPDT[i - 1][j].vnqr);
				Uq[i][j].u[3] = Uq_bar2[i][j].u[3] + Qr.u[3] / 4 * (MPDT[i][j + 1].vnqtheta - 2 * MPDT[i][j].vnqtheta + MPDT[i][j - 1].vnqtheta) + Qz.u[3] / 4 * (MPDT[i + 1][j].vnqtheta - 2 * MPDT[i][j].vnqtheta + MPDT[i - 1][j].vnqtheta);
				Uq[i][j].u[4] = Uq_bar2[i][j].u[4] + Qr.u[4] / 4 * (MPDT[i][j + 1].vnqz - 2 * MPDT[i][j].vnqz + MPDT[i][j - 1].vnqz) + Qz.u[4] / 4 * (MPDT[i + 1][j].vnqz - 2 * MPDT[i][j].vnqz + MPDT[i - 1][j].vnqz);
				Uq[i][j].u[5] = Uq_bar2[i][j].u[5] + Qr.u[5] / 4 * (MPDT[i][j + 1].vpqr - 2 * MPDT[i][j].vpqr + MPDT[i][j - 1].vpqr) + Qz.u[5] / 4 * (MPDT[i + 1][j].vpqr - 2 * MPDT[i][j].vpqr + MPDT[i - 1][j].vpqr);
				Uq[i][j].u[6] = Uq_bar2[i][j].u[6] + Qr.u[6] / 4 * (MPDT[i][j + 1].vpqtheta - 2 * MPDT[i][j].vpqtheta + MPDT[i][j - 1].vpqtheta) + Qz.u[6] / 4 * (MPDT[i + 1][j].vpqtheta - 2 * MPDT[i][j].vpqtheta + MPDT[i - 1][j].vpqtheta);
				Uq[i][j].u[7] = Uq_bar2[i][j].u[7] + Qr.u[7] / 4 * (MPDT[i][j + 1].vpqz - 2 * MPDT[i][j].vpqz + MPDT[i][j - 1].vpqz) + Qz.u[7] / 4 * (MPDT[i + 1][j].vpqz - 2 * MPDT[i][j].vpqz + MPDT[i - 1][j].vpqz);


#ifdef FLUID_DEBUG
				if (i == row && j == col)
				{
					//cout << "final U[i][j].u1 = " << U[i][j].u1 << endl;
					printf("final U[i][j].u1 = %.5e\n", U[i][j].u[0]);
					printf("U_bar2[i][j + 1].u1 = %.5e\n", U_bar2[i][j + 1].u[0]);
					printf("U_bar2[i][j - 1].u1 = %.5e\n", U_bar2[i][j - 1].u[0]);
					printf("U_bar2[i + 1][j].u1 = %.5e\n", U_bar2[i + 1][j].u[0]);
					printf("U_bar2[i - 1][j].u1 = %.5e\n", U_bar2[i - 1][j].u[0]);

					printf("final U[i][j].u3 = %.5e\n", U[i][j].u[2]);
					printf("U_bar2[i][j + 1].u3 = %.5e\n", U_bar2[i][j + 1].u[2]);
					printf("U_bar2[i][j - 1].u3 = %.5e\n", U_bar2[i][j - 1].u[2]);
					printf("U_bar2[i + 1][j].u3 = %.5e\n", U_bar2[i + 1][j].u[2]);
					printf("U_bar2[i - 1][j].u3 = %.5e\n", U_bar2[i - 1][j].u[2]);

				}
#endif// FLUID_DEBUG	

			}
			else if (btype[i][j] == DOWN && ptype[i][j] == CYLINDRICAL_AXIS)
			{
				Qz = arti_q_vis(MPDT[i + 1][j], MPDT[i][j], MPDT[i - 1][j]);

				Uq[i][j].u[0] = Uq_bar2[i][j].u[0] + Qz.u[0] / 4 * (MPDT[i + 1][j].neq - 2 * MPDT[i][j].neq + MPDT[i - 1][j].neq);
				Uq[i][j].u[1] = Uq_bar2[i][j].u[1] + Qz.u[1] / 4 * (MPDT[i + 1][j].peq - 2 * MPDT[i][j].peq + MPDT[i - 1][j].peq);
				Uq[i][j].u[2] = Uq_bar2[i][j].u[2] + Qz.u[2] / 4 * (MPDT[i + 1][j].vnqr - 2 * MPDT[i][j].vnqr + MPDT[i - 1][j].vnqr);
				Uq[i][j].u[3] = Uq_bar2[i][j].u[3] + Qz.u[3] / 4 * (MPDT[i + 1][j].vnqtheta - 2 * MPDT[i][j].vnqtheta + MPDT[i - 1][j].vnqtheta);
				Uq[i][j].u[4] = Uq_bar2[i][j].u[4] + Qz.u[4] / 4 * (MPDT[i + 1][j].vnqz - 2 * MPDT[i][j].vnqz + MPDT[i - 1][j].vnqz);
				Uq[i][j].u[5] = Uq_bar2[i][j].u[5] + Qz.u[5] / 4 * (MPDT[i + 1][j].vpqr - 2 * MPDT[i][j].vpqr + MPDT[i - 1][j].vpqr);
				Uq[i][j].u[6] = Uq_bar2[i][j].u[6] + Qz.u[6] / 4 * (MPDT[i + 1][j].vpqtheta - 2 * MPDT[i][j].vpqtheta + MPDT[i - 1][j].vpqtheta);
				Uq[i][j].u[7] = Uq_bar2[i][j].u[7] + Qz.u[7] / 4 * (MPDT[i + 1][j].vpqz - 2 * MPDT[i][j].vpqz + MPDT[i - 1][j].vpqz);
			}
			else //if (btype[i][j] == 0)
			{
				for (k = 0; k < 13; k++)
				{
					Uq[i][j].u[k] = Uq_bar2[i][j].u[k];
				}
			}
		}
	}

	//更新边界信息



	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			//if (U[i][j].u[0] > 0)
			//{
			//	MPDT[i][j].ne = U[i][j].u[0];
			//}
			//else
			//{
			//	MPDT[i][j].ne = 0;
			//}

			//if (U[i][j].u2 > 0)
			//{
			//	MPDT[i][j].ni = U[i][j].u2;
			//}
			//else
			//{
			//	MPDT[i][j].ni = 0;
			//}

			MPDT[i][j].neq = Uq[i][j].u[0];
			MPDT[i][j].peq = Uq[i][j].u[1];
			MPDT[i][j].vnqr = Uq[i][j].u[2];
			MPDT[i][j].vnqtheta = Uq[i][j].u[3];
			MPDT[i][j].vnqz = Uq[i][j].u[4];


			MPDT[i][j].vpqr = Uq[i][j].u[5];
			MPDT[i][j].vpqtheta = Uq[i][j].u[6];
			MPDT[i][j].vpqz = Uq[i][j].u[7];


		}
	}


	return;

}

struct _U cal_qu(int i, int j)
{
	struct _U uij;
	uij.r = j;
	uij.z = i;


	uij.u[0] = MPDT[i][j].neq;
	uij.u[1] = MPDT[i][j].peq;

	uij.u[2] = MPDT[i][j].vnqr;
	uij.u[3] = MPDT[i][j].vnqtheta;
	uij.u[4] = MPDT[i][j].vnqz;
	uij.u[5] = MPDT[i][j].vpqr;
	uij.u[6] = MPDT[i][j].vpqtheta;
	uij.u[7] = MPDT[i][j].vpqz;

	for (int k = 8; k < 13; k++)
	{
		uij.u[k] = 0;
	}
	return uij;
}

struct _F cal_fqr(struct _U uij)
{
	struct _F fij;

	double neq = uij.u[0];
	double peq = uij.u[1];

	double vnqr = 0;
	double vnqtheta = 0;
	double vnqz = 0;
	double vpqr = 0;
	double vpqtheta = 0;
	double vpqz = 0;

	vnqr = uij.u[2];
	vnqtheta = uij.u[3];
	vnqz = uij.u[4];
	vpqr = uij.u[5];
	vpqtheta = uij.u[6];
	vpqz = uij.u[7];


	fij.r = uij.r;
	fij.z = uij.z;
	fij.f[0] = neq * vnqr;
	fij.f[1] = peq * vpqr;
	for (int k = 2; k < 13; k++)
	{
		fij.f[k] = 0;
	}
	return fij;
}

struct _F cal_fqz(struct _U uij)
{
	struct _F fij;

	double neq = uij.u[0];
	double peq = uij.u[1];

	double vnqr = 0;
	double vnqtheta = 0;
	double vnqz = 0;
	double vpqr = 0;
	double vpqtheta = 0;
	double vpqz = 0;

	vnqr = uij.u[2];
	vnqtheta = uij.u[3];
	vnqz = uij.u[4];
	vpqr = uij.u[5];
	vpqtheta = uij.u[6];
	vpqz = uij.u[7];


	fij.r = uij.r;
	fij.z = uij.z;



	fij.f[0] = neq * vnqz;
	fij.f[1] = peq * vpqz;
	for (int k = 2; k < 13; k++)
	{
		fij.f[k] = 0;
	}
	return fij;
}


struct _F cal_qs(struct _U uij)
{
	struct _F fij;

	double neq = uij.u[0];
	double peq = uij.u[1];

	double vnqr = 0;
	double vnqtheta = 0;
	double vnqz = 0;
	double vpqr = 0;
	double vpqtheta = 0;
	double vpqz = 0;

	vnqr = uij.u[2];
	vnqtheta = uij.u[3];
	vnqz = uij.u[4];
	vpqr = uij.u[5];
	vpqtheta = uij.u[6];
	vpqz = uij.u[7];



	fij.r = uij.r;
	fij.z = uij.z;
	double r1 = 0;
	if (uij.r == 0)
	{
		r1 = -2 / dr;
	}
	else
	{
		r1 = -1 / (uij.r * dr);
	}
	fij.f[0] = r1 * neq * vnqr;
	fij.f[1] = r1 * peq * vpqr;

	for (int k = 2; k < 13; k++)
	{
		fij.f[k] = 0;
	}
	return fij;
}

struct _U arti_q_vis(struct  node np, struct  node n, struct  node nn)
{
	struct _U Q;
	double eta = 0.05;

	for (int k = 0; k < 13; k++)
	{
		Q.u[k] = 0;
	}

	if (abs(np.neq - n.neq) + abs(n.neq - nn.neq) != 0)
	{
		Q.u[0] = eta * abs(abs(np.neq - n.neq) - abs(n.neq - nn.neq)) / abs(abs(np.neq - n.neq) + abs(n.neq - nn.neq));
	}

	if (abs(np.peq - n.peq) + abs(n.peq - nn.peq) != 0)
	{
		Q.u[1] = eta * abs(abs(np.peq - n.peq) - abs(n.peq - nn.peq)) / abs(abs(np.peq - n.peq) + abs(n.peq - nn.peq));
	}

	Q.u[2] = Q.u[0];
	Q.u[3] = Q.u[0];
	Q.u[4] = Q.u[0];
	Q.u[11] = Q.u[0];

	Q.u[5] = Q.u[1];
	Q.u[6] = Q.u[1];
	Q.u[7] = Q.u[1];
	Q.u[12] = Q.u[1];

	return Q;
}


void smooth_ne(int i, int j)
{
	int is, js, je, ie;

	for (i = 0; i < nz; i++)
	{
		for (j = 0; j < nr; j++)
		{
			if (MPDT[i][j].ne < 0)
			{
				if (i == 0)
				{
					is = i;
					ie = i + 2;
				}
				else if (i == nz - 1)
				{
					ie = i;
					is = i - 2;
				}
				else
				{
					ie = i + 1;
					is = i - 1;
				}

				if (j == 0)
				{
					js = j;
					je = j + 2;
				}
				else if (i == nz - 1)
				{
					js = j - 2;
					je = j;
				}
				else
				{
					js = j - 1;
					je = j + 1;
				}
				double tsum = 0;
				double cn = 0;
				for (int id = is; id <= ie; id++)
				{
					for (int jd = js; jd <= je; jd++)
					{
						if (MPDT[id][jd].ne > 0)
						{
							tsum += MPDT[id][jd].ne;
							cn++;
						}
					}
				}

				MPDT[i][j].ne = tsum / cn;
			}
		}
	}
}

void smooth_pe(int i, int j)
{
	int is, js, je, ie;

	if (MPDT[i][j].pe < 0)
	{
		if (i == 0)
		{
			is = i;
			ie = i + 2;
		}
		else if (i == nz - 1)
		{
			ie = i;
			is = i - 2;
		}
		else
		{
			ie = i + 1;
			is = i - 1;
		}

		if (j == 0)
		{
			js = j;
			je = j + 2;
		}
		else if (i == nz - 1)
		{
			js = j - 2;
			je = j;
		}
		else
		{
			js = j - 1;
			je = j + 1;
		}
		double tsum = 0;
		double cn = 0;
		for (int id = is; id <= ie; id++)
		{
			for (int jd = js; jd <= je; jd++)
			{
				if (MPDT[id][jd].pe > 0)
				{
					tsum += MPDT[id][jd].pe;
					cn++;
				}
			}
		}

		MPDT[i][j].pe = tsum / cn;
	}
}