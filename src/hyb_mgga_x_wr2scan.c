/*
 Copyright (C) 2020 Susi Lehtola

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "util.h"

#define XC_HYB_MGGA_X_WR2SCAN    662 /* Hybrid version of re-regularized SCAN exchange */

typedef struct{
  double c1, c2, d, k1;
  double eta, dp2;
  double alpha;
  double beta;
  double omega;
} mgga_x_wr2scan_params;

#define N_PAR 9
static const char *names[N_PAR] = {"_c1", "_c2", "_d", "_k1", "_eta", "_dp2", "_alpha", "_beta", "_omega"};
static const char *desc[N_PAR] = {"c1 parameter", "c2 parameter", "d parameter",
                                  "k1 parameter", "eta parameter", "dp2 parameter",
                                  "fraction of full Hartree-Fock exchange",
                                  "fraction of short-range exact exchange",
                                  "range separation parameter"};

static const double par_wr2scan[N_PAR] = {0.667, 0.8, 1.24, 0.065, 0.001, 0.361, 1.0, -1.0, 0.3};

static void
hyb_mgga_x_wr2scan_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_wr2scan_params));
  xc_hyb_init_cam(p, 0.0, 0.0, 0.0);
}

#include "maple2c/mgga_exc/hyb_mgga_x_wr2scan.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_hyb_mgga_x_wr2scan = {
  XC_HYB_MGGA_X_WR2SCAN,
  XC_EXCHANGE,
  "Range-separated re-regularized SCAN exchange by Wittmann et al",
  XC_FAMILY_HYB_MGGA,
  {&xc_ref_Wittmann2023_224103, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_HYB_CAM | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-11,
  {N_PAR, names, desc, par_wr2scan, set_ext_params_cpy_cam},
  hyb_mgga_x_wr2scan_init, NULL,
  NULL, NULL, &work_mgga
};
