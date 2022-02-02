#pragma once

#include "MHD.h"


void electron_flow();
void electron_flow_v2();
void ion_flow();

struct _U cal_u(int i, int j);
struct _F cal_fr(struct _U uij);
struct _F cal_fz(struct _U uij);
struct _F cal_s(struct _U uij);
struct _U arti_vis(struct  node np, struct  node n, struct  node nn);

void Q_fluid();
struct _U cal_qu(int i, int j);
struct _F cal_fqr(struct _U uij);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
struct _F cal_fqz(struct _U uij);
struct _F cal_qs(struct _U uij);

struct _U arti_q_vis(struct  node np, struct  node n, struct  node nn);

void smooth_ne(int i, int j);
void smooth_pe(int i, int j);