#include "MHD.h"

/*************************************************

Copyright (c) ,2015-2016,  
File name:
Author:   

Date:2010-08-25

Description:描述主要实现的功能
Fucntion list:...

Histiory : //修改历史记录列表，每条修改记录修改日期，修改者，以及修改内容

1.Data:....
  Author:...
  Midification:...

2.


**************************************************/



/*************************************************

Function:                  // 函数名称

Description:               // hanshu ogngnegn 

Calls:                     // 无

Input:                     // 输入参数说明，包括每个参数的作

						   // 用、取值说明及参数间关系。

Output:                    // 对输出参数的说明。

Return:                    // 函数返回值的说明

Others:                   // 其它说明

*************************************************/

double max_q_speed = 1e-5;

void move()
{
	double q_half;
	double hrho, htheta, hz, h2;
	double srho, stheta, sz;
	double urho, utheta, uz;
	double urho_half, utheta_half, uz_half;


	for (int i = 0; i < nz; i++)
	{
		for (int j = 0; j < nr; j++)
		{
			//if (i == 390 && j == 0)
			//{
			//	printf("ptype[i][j] = %d\n", ptype[i][j]);
			//}
			if ((ptype[i][j] & VACCUM_BOUNDARY) != 0) continue;
			if (MPDT[i][j].ne == 0 || MPDT[i][j].ni == 0) continue;
			if ((ptype[i][j] & ANODE_BOUNDARY) != 0 || (ptype[i][j] & CATHODE_BOUNDARY) != 0 ) continue;
			q_half = -dt * QE / ME / 2;
			hrho = q_half * app_Br[i][j];
			htheta = 0;
			hz = q_half * app_Bz[i][j];
			h2 = sqr(hrho) + sqr(hz);
			
			srho = 2 * hrho / (1 + h2);
			stheta = 2 * htheta / (1 + h2); 
			sz = 2 * hz / (1 + h2);
			urho = MPDT[i][j].ver + q_half * Er[i][j];
			utheta = MPDT[i][j].vetheta;
			uz = MPDT[i][j].vez + q_half * Ez[i][j];
			urho_half = urho + (utheta * sz - uz * stheta) + ((uz * hrho - urho * hz) * sz - (urho * htheta - utheta * hrho) * stheta);
			utheta_half = utheta + (uz * srho - urho * sz) + ((urho * htheta - utheta * hrho) * srho - (utheta * hz - uz * htheta) * sz);
			uz_half = uz + (urho * stheta - utheta * srho) + ((utheta * hz - uz * htheta) * stheta - (uz * hrho - urho * hz) * srho);
			MPDT[i][j].ver = urho_half + q_half * Er[i][j];
			MPDT[i][j].vetheta = utheta_half;
			MPDT[i][j].vez = uz_half + q_half * Ez[i][j];
		
			//MPDT[i][j].ver = -QE * Er[i][j] / ME * dt;
			//MPDT[i][j].vez = -QE * Ez[i][j] / ME * dt;

			//电磁场加速之后电子流体动能密度
			double pre_ee   = MPDT[i][j].pe /(MPDT[i][j].ne * (gamma - 1) ) + 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);

			//if (MPDT[i][j].ver > 1.8e6)
			//{
			//	MPDT[i][j].ver = 1.8e6;
			//}
			//else if (MPDT[i][j].ver < -1.8e6)
			//{
			//	MPDT[i][j].ver = -1.8e6;
			//}

			//if (abs(MPDT[i][j].vetheta) > 1.8e6)
			//{
			//	MPDT[i][j].vetheta = 1.8e6;
			//}
			//else if (MPDT[i][j].vetheta < -1.8e6)
			//{
			//	MPDT[i][j].vetheta = -1.8e6;
			//}

			//if (MPDT[i][j].vez > 1.8e6)
			//{
			//	MPDT[i][j].vez = 1.8e6;

			//}
			//else if (MPDT[i][j].vez < -1.8e6)
			//{
			//	MPDT[i][j].vez = -1.8e6;
			//}


			q_half = dt * QE / MI / 2;
			hrho = q_half * app_Br[i][j];
			htheta = 0;
			hz = q_half * app_Bz[i][j];
			h2 = sqr(hrho) + sqr(hz);

			srho = 2 * hrho / (1 + h2);
			stheta = 2 * htheta / (1 + h2);
			sz = 2 * hz / (1 + h2);
			urho = MPDT[i][j].vir + q_half * Er[i][j];
			utheta = MPDT[i][j].vitheta;
			uz = MPDT[i][j].viz + q_half * Ez[i][j];
			urho_half = urho + (utheta * sz - uz * stheta) + ((uz * hrho - urho * hz) * sz - (urho * htheta - utheta * hrho) * stheta);
			utheta_half = utheta + (uz * srho - urho * sz) + ((urho * htheta - utheta * hrho) * srho - (utheta * hz - uz * htheta) * sz);
			uz_half = uz + (urho * stheta - utheta * srho) + ((utheta * hz - uz * htheta) * stheta - (uz * hrho - urho * hz) * srho);
			MPDT[i][j].vir = urho_half + q_half * Er[i][j];
			MPDT[i][j].vitheta = utheta_half;
			MPDT[i][j].viz = uz_half + q_half * Ez[i][j];

			coulomb_collision(i, j);
			ionization_collisions(i, j);
		//	double pre_ei = 0.5 * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);



		//	double Q_ei = 1.4e-22;//电子离子碰撞截面
		//	MPDT[i][j].tau_ei = pow(pre_ee,1.5) * sqrt(ME) / (11.313708 * PI * MPDT[i][j].ni * pow(QE, 4) * 10);
		//	//double sigma_q = PI * pow(QE, 4) / (128 * sqr(EPS_0) * sqr(ME) * pow((MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz), 4));
		//	MPDT[i][j].sigma_Q = PI * pow(QE, 4) / (128 * sqr(EPS_0) * sqr(ME) * pow((MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz), 4));
		//	MPDT[i][j].mu_ie = MPDT[i][j].ni * Q_ei * sqrt((12 * MPDT[i][j].ee) / (PI * ME));

		//	//if(abs(MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz) < 1e-6)
		//	//{ 
		//	//	MPDT[i][j].mu_ie = 0;
		//	//}
		//	//else
		//	//{
		//	//	MPDT[i][j].mu_ie = (MPDT[i][j].ne * PI * pow(QE, 4)) / (128 * sqr(EPS_0) * sqr(ME) * pow((MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz), 3));
		//	//}
		//	double m_eir = 0;        //电子向离子转移的径向（r方向）动量
		//	double m_eitheta = 0;    //电子向离子转移的角向（theta方向）动量
		//	double m_eiz = 0;        //电子向离子转移的轴向（z方向）动量

		//	MPDT[i][j].delta_ei = dt / MPDT[i][j].tau_ei * (pre_ee - pre_ei);     //电子向离子转移的能量
		//	
		//	//double eta = MPDT[i][j].mu_ie * MPDT[i][j].ne * ME;
		//	//m_eir      = eta * (MPDT[i][j].ver - MPDT[i][j].vir);
		//	//m_eitheta  = eta * (MPDT[i][j].vetheta - MPDT[i][j].vitheta);
		//	//m_eiz      = eta * (MPDT[i][j].vez - MPDT[i][j].viz);

		//	////delta_ei = 3 * MPDT[i][j].ne * ME * MPDT[i][j].mu_ie * (MPDT[i][j].ee - MPDT[i][j].ei);

		//	//double pre_eng = MPDT[i][j].ei + MPDT[i][j].ee;

		//	//if (MPDT[i][j].ne != 0 && MPDT[i][j].ni != 0)
		//	//{
		//	//	MPDT[i][j].ver -= m_eir / (MPDT[i][j].ne * ME);
		//	//	MPDT[i][j].vir += m_eir / (MPDT[i][j].ni * MI);
		//	//	MPDT[i][j].vetheta -= m_eitheta / (MPDT[i][j].ne * ME);
		//	//	MPDT[i][j].vitheta += m_eitheta / (MPDT[i][j].ni * MI);
		//	//	MPDT[i][j].vez -= m_eiz / (MPDT[i][j].ne * ME);
		//	//	MPDT[i][j].viz += m_eiz / (MPDT[i][j].ni * MI);
		//	//}
		//	tur = sqr(MPDT[i][j].ver - MPDT[i][j].vir);
		//	tutheta = sqr(MPDT[i][j].vetheta - MPDT[i][j].vitheta);
		//	tuz = sqr(MPDT[i][j].vez - MPDT[i][j].viz);

		//	tu2 = tuz + tutheta + tuz;

		//	MPDT[i][j].vir += sqrt(2 * tur / tu2 * MPDT[i][j].delta_ei / (MPDT[i][j].ni * MI));
		//	MPDT[i][j].vitheta += sqrt(2 * tutheta / tu2 * MPDT[i][j].delta_ei / (MPDT[i][j].ni * MI));
		//	MPDT[i][j].viz += sqrt(2 * tuz / tu2 * MPDT[i][j].delta_ei / (MPDT[i][j].ni * MI));

		//	MPDT[i][j].ver -= sqrt(2 * tur / tu2 * 11 * MPDT[i][j].delta_ei / (MPDT[i][j].ne * ME));
		//	MPDT[i][j].vetheta -= sqrt(2 * tutheta / tu2 * 11 * MPDT[i][j].delta_ei / (MPDT[i][j].ne * ME));
		//	MPDT[i][j].vez -= sqrt(2 * tuz / tu2 * 11 * MPDT[i][j].delta_ei / (MPDT[i][j].ne * ME));

		//	//MPDT[i][j].ver = MPDT[i][j].vir;
		//	//MPDT[i][j].vetheta = MPDT[i][j].vitheta;
		//	//MPDT[i][j].vez = MPDT[i][j].viz;

		//	//离子速度和磁场夹角
	MPDT[i][j].angle_b_vi = magnetic_vec_angle(app_Br[i][j], app_Bz[i][j], MPDT[i][j].vir, MPDT[i][j].viz);

		//	////碰撞之后电子流体能量
		//	//double alter_ee = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);

		//	//double pre_ei = MPDT[i][j].ei;
		//	//double ev = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez) ;
		//	MPDT[i][j].ei = 0.5 * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz) ;

		//	//double dee = pre_ee - alter_ee;
		//	//double dei = MPDT[i][j].ei - pre_ei;


		//	MPDT[i][j].pe += MPDT[i][j].ne * 10 * MPDT[i][j].delta_ei * (gamma - 1);
		//	MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne) + 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);




		//	//电离部分待完成，电子能量大于电离值增加电子离子密度

		//	if (MPDT[i][j].ee > 15.6 * QE)
		//	{
		//		if (MPDT[i][j].ne < 1e21)
		//		{
		//			MPDT[i][j].ni += MPDT[i][j].ne;
		//			MPDT[i][j].ne *= 2;
		//		}

		//		double tep = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne);
		//		double teu = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);

		//		double tep1 = tep / (tep + teu) * 15.6 * QE;
		//		double teu2 = teu / (tep + teu) * 15.6 * QE;


		//		tur = sqr(MPDT[i][j].ver);
		//		tutheta = sqr(MPDT[i][j].vetheta);
		//		tuz = sqr(MPDT[i][j].vez);

		//		tu2 = tuz + tutheta + tuz;

		//		MPDT[i][j].ver -= sqrt(2 * tur / tu2 * teu2 /(MPDT[i][j].ne * ME));
		//		MPDT[i][j].vetheta -= sqrt(2 * tutheta / tu2 * teu2 / (MPDT[i][j].ne * ME));
		//		MPDT[i][j].vez -= sqrt(2 * tuz / tu2 * teu2 / (MPDT[i][j].ne * ME));

		//		MPDT[i][j].pe -= MPDT[i][j].ne * tep1 * (gamma - 1);
		//		MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne) + 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
		//	}
		}
	}
}

void move_q()
{
	double q_half;
	double hrho, htheta, hz, h2;
	double srho, stheta, sz;
	double urho, utheta, uz;
	double urho_half, utheta_half, uz_half;
	for (int i = 0; i < nz; i++)
	{
		for (int j = 0; j < nr; j++)
		{
			if ((ptype[i][j] & VACCUM_BOUNDARY) != 0) continue;
			if (MPDT[i][j].neq != 0 ) 
			{
				//判断电子电荷和离子电荷是否分离
				if (is_electron_ion_separation(MPDT[i][j].angle_b_vi) == true)
				{
					//计算多余电子运动
					q_half = -dt * QE / ME / 2;
					hrho = q_half * app_Br[i][j];
					htheta = 0;
					hz = q_half * app_Bz[i][j];
					h2 = sqr(hrho) + sqr(hz);

					srho = 2 * hrho / (1 + h2);
					stheta = 2 * htheta / (1 + h2);
					sz = 2 * hz / (1 + h2);
					urho = MPDT[i][j].vnqr + q_half * Er[i][j];
					utheta = MPDT[i][j].vnqtheta;
					uz = MPDT[i][j].vnqz + q_half * Ez[i][j];
					urho_half = urho + (utheta * sz - uz * stheta) + ((uz * hrho - urho * hz) * sz - (urho * htheta - utheta * hrho) * stheta);
					utheta_half = utheta + (uz * srho - urho * sz) + ((urho * htheta - utheta * hrho) * srho - (utheta * hz - uz * htheta) * sz);
					uz_half = uz + (urho * stheta - utheta * srho) + ((utheta * hz - uz * htheta) * stheta - (uz * hrho - urho * hz) * srho);
					MPDT[i][j].vnqr = urho_half + q_half * Er[i][j];
					MPDT[i][j].vnqtheta = utheta_half;
					MPDT[i][j].vnqz = uz_half + q_half * Ez[i][j];

					urho = MPDT[i][j].vnqr;
					utheta = MPDT[i][j].vnqtheta;
					uz = MPDT[i][j].vnqz;


					if (is_large_max_speed(urho, utheta, uz, max_q_speed) == true)
					{
						double tu = sqrt(sqr(urho) + sqr(utheta) + sqr(uz));
						MPDT[i][j].vnqr = urho / tu * max_q_speed;
						MPDT[i][j].vnqtheta = utheta / tu * max_q_speed;
						MPDT[i][j].vnqz = uz / tu * max_q_speed;
					}
				}
				else
				{
					if (is_large_max_speed(MPDT[i][j].vir, MPDT[i][j].vitheta, MPDT[i][j].viz, max_q_speed) == true)
					{
						double tu = sqrt(sqr(MPDT[i][j].vir) + sqr(MPDT[i][j].vitheta) + sqr(MPDT[i][j].viz));
						MPDT[i][j].vnqr = MPDT[i][j].vir / tu * max_q_speed;
						MPDT[i][j].vnqtheta = MPDT[i][j].vitheta / tu * max_q_speed;
						MPDT[i][j].vnqz = MPDT[i][j].viz / tu * max_q_speed;
					}
					else
					{
						MPDT[i][j].vnqr = MPDT[i][j].vir;
						MPDT[i][j].vnqtheta = MPDT[i][j].vitheta;
						MPDT[i][j].vnqz = MPDT[i][j].viz;
					}
				}
			}

			if (MPDT[i][j].peq != 0)
			{
				q_half = dt * QE / ME / 2;
				hrho = q_half * app_Br[i][j];
				htheta = 0;
				hz = q_half * app_Bz[i][j];
				h2 = sqr(hrho) + sqr(hz);

				srho = 2 * hrho / (1 + h2);
				stheta = 2 * htheta / (1 + h2);
				sz = 2 * hz / (1 + h2);
				urho = MPDT[i][j].vpqr + q_half * Er[i][j];
				utheta = MPDT[i][j].vpqtheta;
				uz = MPDT[i][j].vpqz + q_half * Ez[i][j];
				urho_half = urho + (utheta * sz - uz * stheta) + ((uz * hrho - urho * hz) * sz - (urho * htheta - utheta * hrho) * stheta);
				utheta_half = utheta + (uz * srho - urho * sz) + ((urho * htheta - utheta * hrho) * srho - (utheta * hz - uz * htheta) * sz);
				uz_half = uz + (urho * stheta - utheta * srho) + ((utheta * hz - uz * htheta) * stheta - (uz * hrho - urho * hz) * srho);
				MPDT[i][j].vpqr = urho_half + q_half * Er[i][j];
				MPDT[i][j].vpqtheta = utheta_half;
				MPDT[i][j].vpqz = uz_half + q_half * Ez[i][j];


				urho = MPDT[i][j].vpqr;
				utheta = MPDT[i][j].vpqtheta;
				uz = MPDT[i][j].vpqz;

				if (is_large_max_speed(urho, utheta, uz, max_q_speed) == true)
				{
					double tu = sqrt(sqr(urho) + sqr(utheta) + sqr(uz));
					MPDT[i][j].vpqr = urho / tu * max_q_speed;
					MPDT[i][j].vpqtheta = utheta / tu * max_q_speed;
					MPDT[i][j].vpqz = uz / tu * max_q_speed;
				}
			}
			
		}
	}
}

void plasma_para()
{

}


void collision()
{
	
}

double  magnetic_vec_angle(double var,double vaz,double vbr,double vbz)
{
	double tmp = (var * vbr + vaz * vbz) / (sqrt(sqr(var) + sqr(vaz)) * sqrt(sqr(vbr) + sqr(vbz)));
	return acos(tmp) * 180 / PI;
}

bool is_electron_ion_separation(double angle)
{
	return true;
	//return angle > 10;
}

bool is_large_max_speed(double ur, double utheta, double uz, double max_speed)
{
	return sqr(ur) + sqr(utheta) + sqr(uz) > sqr(max_speed);
}

//double ionization_collisions(int i, int j)
//{
//	double rend = 0;
//
//	return 0;
//}

void ionization_collisions(int i, int j)
{
	double tur;
	double tutheta;
	double tuz;

	double tu2;
	double nn = (1e21 - MPDT[i][j].ne - MPDT[i][j].neq + MPDT[i][j].peq);

	if (nn < 0)
	{
		return;
	}

	double nee = MPDT[i][j].ee + 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
	double p = 1 - exp( -nn * sqrt(2 * nee / ME / 3) * 1e-20 * dt);

	if (MPDT[i][j].ee - 15.755 * QE * p > 1e-25)
	{
		//MPDT[i][j].ni += nn * p;
		//MPDT[i][j].ne += nn * p;

		MPDT[i][j].ee -= 15.755 * QE * p;

		tur = sqr(MPDT[i][j].ver);
		tutheta = sqr(MPDT[i][j].vetheta);
		tuz = sqr(MPDT[i][j].vez);

		//tu2 = tuz + tutheta + tuz;

		//MPDT[i][j].ver -= sqrt(2 * tur / tu2 * teu2 / (MPDT[i][j].ne * ME));
		//MPDT[i][j].vetheta -= sqrt(2 * tutheta / tu2 * teu2 / (MPDT[i][j].ne * ME));
		//MPDT[i][j].vez -= sqrt(2 * tuz / tu2 * teu2 / (MPDT[i][j].ne * ME));

		//double add_n = nn * p;

		//MPDT[i][j].ver = MPDT[i][j].ver * MPDT[i][j].ne / (MPDT[i][j].ne + add_n);
		//MPDT[i][j].vetheta = MPDT[i][j].vetheta * MPDT[i][j].ne / (MPDT[i][j].ne + add_n);
		//MPDT[i][j].vez = MPDT[i][j].vez * MPDT[i][j].ne / (MPDT[i][j].ne + add_n);

		//MPDT[i][j].vir = MPDT[i][j].vir * MPDT[i][j].ne / (MPDT[i][j].ni + add_n);
		//MPDT[i][j].vitheta = MPDT[i][j].vitheta * MPDT[i][j].ne / (MPDT[i][j].ni + add_n);
		//MPDT[i][j].viz = MPDT[i][j].viz * MPDT[i][j].ne / (MPDT[i][j].ni + add_n);

		MPDT[i][j].ni += nn * p;
		MPDT[i][j].ne += nn * p;

		//MPDT[i][j].pe -= MPDT[i][j].ne * tep1 * (gamma - 1);
		//MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne) + 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
		//printf("ee pre = %e,after = %e\n", nee, MPDT[i][j].ee);
		MPDT[i][j].ei = 0.5 * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);
	}

	/*if (MPDT[i][j].ee > 15.755 * QE)
	{

		MPDT[i][j].ni += MPDT[i][j].ne;
		MPDT[i][j].ne *= 2;

		double tep = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne);
		double teu = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);

		double tep1 = tep / (tep + teu) * 15.6 * QE;
		double teu2 = teu / (tep + teu) * 15.6 * QE;


		tur = sqr(MPDT[i][j].ver);
		tutheta = sqr(MPDT[i][j].vetheta);
		tuz = sqr(MPDT[i][j].vez);

		tu2 = tuz + tutheta + tuz;

		MPDT[i][j].ver -= sqrt(2 * tur / tu2 * teu2 / (MPDT[i][j].ne * ME));
		MPDT[i][j].vetheta -= sqrt(2 * tutheta / tu2 * teu2 / (MPDT[i][j].ne * ME));
		MPDT[i][j].vez -= sqrt(2 * tuz / tu2 * teu2 / (MPDT[i][j].ne * ME));

		MPDT[i][j].pe -= MPDT[i][j].ne * tep1 * (gamma - 1);
		MPDT[i][j].ee = MPDT[i][j].pe / (gamma - 1) / (MPDT[i][j].ne) + 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
	}*/

	return;
}

/// <summary>
/// 库仑碰撞，进行电子离子动量和动能转移
/// </summary>
/// <param name="i"></param>
/// <param name="j"></param>
/// <returns></returns>
void coulomb_collision(int i, int j)
{
	double ur, utheta, uz;
	double ur1, utheta1, uz1;
	double pre_ee, pre_ei;
	double ep;
	double ne = MPDT[i][j].ne + MPDT[i][j].neq  - MPDT[i][j].peq;

	pre_ee = 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
	pre_ei = 0.5 * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);

	ur = (MPDT[i][j].ver * ne * ME + MPDT[i][j].vir * MPDT[i][j].ni * MI) / (ne * ME + MPDT[i][j].ni * MI);
	utheta = (MPDT[i][j].vetheta * ne * ME + MPDT[i][j].vitheta * MPDT[i][j].ni * MI) / (ne * ME + MPDT[i][j].ni * MI);
	uz = (MPDT[i][j].vez * ne * ME + MPDT[i][j].viz * MPDT[i][j].ni * MI) / (ne * ME + MPDT[i][j].ni * MI);

	ep = pre_ee + pre_ei - 0.5 * (ME + MI) * (sqr(ur) + sqr(utheta) + sqr(uz));
	//MPDT[i][j].pe = ne * ep * (gamma - 1);

	ur1 = sqrt((sqr(ur) + sqr(utheta) + sqr(uz)) / 3.0);

	MPDT[i][j].vir = ur + 0.01 * (ur1 - ur);
	MPDT[i][j].vitheta = utheta + 0.01 * (ur1 - utheta);
	MPDT[i][j].viz = uz + 0.01 * (ur1 - uz);

	MPDT[i][j].ver = MPDT[i][j].vir;
	MPDT[i][j].vetheta = MPDT[i][j].vitheta;
	MPDT[i][j].vez = MPDT[i][j].viz;

	//MPDT[i][j].vir = ur;
	//MPDT[i][j].vitheta = utheta;
	//MPDT[i][j].viz = uz;

	//MPDT[i][j].ver = ur;
	//MPDT[i][j].vetheta = utheta;
	//MPDT[i][j].vez = uz;
	
	MPDT[i][j].ee += ep;

	if (MPDT[i][j].ee > 180.7 * QE)
	{
		MPDT[i][j].ee = 180.7 * QE;
	}

	//MPDT[i][j].ee = ep + 0.5 * ME * (MPDT[i][j].ver * MPDT[i][j].ver + MPDT[i][j].vetheta * MPDT[i][j].vetheta + MPDT[i][j].vez * MPDT[i][j].vez);
	MPDT[i][j].ei = 0.5 * MI * (MPDT[i][j].vir * MPDT[i][j].vir + MPDT[i][j].vitheta * MPDT[i][j].vitheta + MPDT[i][j].viz * MPDT[i][j].viz);
	return;
}

void recombination_collision()
{
	
}

void radiation()
{

}

void current_caulate()
{
	int start_r = 0, end_r = 0, pos_z = 0;
	int start_z = 0, end_z = 0, pos_r = 0;
	register int i, j;

	pos_r = 34;
	start_z = 0;
	end_z = 382;

	pos_z = 382;
	start_r = 0;
	end_r = 34;

	current_I = 0;

	j = pos_r;

	for(i = start_z; i < end_z; i++)
	{
		current_I += (MPDT[i][j].peq * MPDT[i][j].vpqr - MPDT[i][j].neq * MPDT[i][j].vnqr) * QE  * dr * dtheta * j * dz /dt;
		//current_I += (MPDT[i][j].peq * MPDT[i][j].vpqz - MPDT[i][j].neq * MPDT[i][j].vnqz) * QE ;
	}
	i = pos_z;

	for (j = start_r; j < end_r; j++)
	{
		//current_I += (MPDT[i][j].peq * MPDT[i][j].vpqr - MPDT[i][j].neq * MPDT[i][j].vnqr) * QE 
		if (j != 0) 
		{
			current_I += (MPDT[i][j].peq * MPDT[i][j].vpqz - MPDT[i][j].neq * MPDT[i][j].vnqz) * QE  * sqr(dr) * dtheta * (2 * j - 1) / dt;
		}
		
	}

	printf("current = %lf\n", current_I);
}


void current_control()
{
	pid.err = pid.set_current - (-current_I);
	pid.integral += pid.err;
	pid.ne_density = pid.Kp * pid.err + pid.Ki * pid.integral + pid.Kd * (pid.err - pid.err_last);

	pid.err_last = pid.err;

	inter_e_den += pid.ne_density;
	if (inter_e_den < 0)
	{
		inter_e_den = 0;
	}
	out_e_den = inter_e_den * cathod_cell / anode_cell;

}