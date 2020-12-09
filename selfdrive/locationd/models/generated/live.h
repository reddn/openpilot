/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_9092064331387050047);
void inv_err_fun(double *nom_x, double *true_x, double *out_5072860975779674511);
void H_mod_fun(double *state, double *out_1387527780010964011);
void f_fun(double *state, double dt, double *out_472167469385665616);
void F_fun(double *state, double dt, double *out_843513716817994442);
void h_3(double *state, double *unused, double *out_4650146831832270943);
void H_3(double *state, double *unused, double *out_726371485440298204);
void h_4(double *state, double *unused, double *out_5900747621767450618);
void H_4(double *state, double *unused, double *out_720594762515236831);
void h_9(double *state, double *unused, double *out_6531791413608352675);
void H_9(double *state, double *unused, double *out_8357614478688601990);
void h_10(double *state, double *unused, double *out_2636647597143245490);
void H_10(double *state, double *unused, double *out_4171743780302226264);
void h_12(double *state, double *unused, double *out_6755562667885504704);
void H_12(double *state, double *unused, double *out_2827150274254446400);
void h_31(double *state, double *unused, double *out_2817406634413468316);
void H_31(double *state, double *unused, double *out_152060499111026579);
void h_32(double *state, double *unused, double *out_1225412177397073175);
void H_32(double *state, double *unused, double *out_1798441222536725068);
void h_13(double *state, double *unused, double *out_7788748001042336296);
void H_13(double *state, double *unused, double *out_9177661506835821808);
void h_14(double *state, double *unused, double *out_6531791413608352675);
void H_14(double *state, double *unused, double *out_8357614478688601990);
void h_19(double *state, double *unused, double *out_4261745000301853400);
void H_19(double *state, double *unused, double *out_3924934627586170809);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);