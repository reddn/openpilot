
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5091353622417362108) {
   out_5091353622417362108[0] = delta_x[0] + nom_x[0];
   out_5091353622417362108[1] = delta_x[1] + nom_x[1];
   out_5091353622417362108[2] = delta_x[2] + nom_x[2];
   out_5091353622417362108[3] = delta_x[3] + nom_x[3];
   out_5091353622417362108[4] = delta_x[4] + nom_x[4];
   out_5091353622417362108[5] = delta_x[5] + nom_x[5];
   out_5091353622417362108[6] = delta_x[6] + nom_x[6];
   out_5091353622417362108[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3347397060309790859) {
   out_3347397060309790859[0] = -nom_x[0] + true_x[0];
   out_3347397060309790859[1] = -nom_x[1] + true_x[1];
   out_3347397060309790859[2] = -nom_x[2] + true_x[2];
   out_3347397060309790859[3] = -nom_x[3] + true_x[3];
   out_3347397060309790859[4] = -nom_x[4] + true_x[4];
   out_3347397060309790859[5] = -nom_x[5] + true_x[5];
   out_3347397060309790859[6] = -nom_x[6] + true_x[6];
   out_3347397060309790859[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_8038179413808804489) {
   out_8038179413808804489[0] = 1.0;
   out_8038179413808804489[1] = 0.0;
   out_8038179413808804489[2] = 0.0;
   out_8038179413808804489[3] = 0.0;
   out_8038179413808804489[4] = 0.0;
   out_8038179413808804489[5] = 0.0;
   out_8038179413808804489[6] = 0.0;
   out_8038179413808804489[7] = 0.0;
   out_8038179413808804489[8] = 0.0;
   out_8038179413808804489[9] = 1.0;
   out_8038179413808804489[10] = 0.0;
   out_8038179413808804489[11] = 0.0;
   out_8038179413808804489[12] = 0.0;
   out_8038179413808804489[13] = 0.0;
   out_8038179413808804489[14] = 0.0;
   out_8038179413808804489[15] = 0.0;
   out_8038179413808804489[16] = 0.0;
   out_8038179413808804489[17] = 0.0;
   out_8038179413808804489[18] = 1.0;
   out_8038179413808804489[19] = 0.0;
   out_8038179413808804489[20] = 0.0;
   out_8038179413808804489[21] = 0.0;
   out_8038179413808804489[22] = 0.0;
   out_8038179413808804489[23] = 0.0;
   out_8038179413808804489[24] = 0.0;
   out_8038179413808804489[25] = 0.0;
   out_8038179413808804489[26] = 0.0;
   out_8038179413808804489[27] = 1.0;
   out_8038179413808804489[28] = 0.0;
   out_8038179413808804489[29] = 0.0;
   out_8038179413808804489[30] = 0.0;
   out_8038179413808804489[31] = 0.0;
   out_8038179413808804489[32] = 0.0;
   out_8038179413808804489[33] = 0.0;
   out_8038179413808804489[34] = 0.0;
   out_8038179413808804489[35] = 0.0;
   out_8038179413808804489[36] = 1.0;
   out_8038179413808804489[37] = 0.0;
   out_8038179413808804489[38] = 0.0;
   out_8038179413808804489[39] = 0.0;
   out_8038179413808804489[40] = 0.0;
   out_8038179413808804489[41] = 0.0;
   out_8038179413808804489[42] = 0.0;
   out_8038179413808804489[43] = 0.0;
   out_8038179413808804489[44] = 0.0;
   out_8038179413808804489[45] = 1.0;
   out_8038179413808804489[46] = 0.0;
   out_8038179413808804489[47] = 0.0;
   out_8038179413808804489[48] = 0.0;
   out_8038179413808804489[49] = 0.0;
   out_8038179413808804489[50] = 0.0;
   out_8038179413808804489[51] = 0.0;
   out_8038179413808804489[52] = 0.0;
   out_8038179413808804489[53] = 0.0;
   out_8038179413808804489[54] = 1.0;
   out_8038179413808804489[55] = 0.0;
   out_8038179413808804489[56] = 0.0;
   out_8038179413808804489[57] = 0.0;
   out_8038179413808804489[58] = 0.0;
   out_8038179413808804489[59] = 0.0;
   out_8038179413808804489[60] = 0.0;
   out_8038179413808804489[61] = 0.0;
   out_8038179413808804489[62] = 0.0;
   out_8038179413808804489[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_6958450483037231644) {
   out_6958450483037231644[0] = state[0];
   out_6958450483037231644[1] = state[1];
   out_6958450483037231644[2] = state[2];
   out_6958450483037231644[3] = state[3];
   out_6958450483037231644[4] = state[4];
   out_6958450483037231644[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6958450483037231644[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6958450483037231644[7] = state[7];
}
void F_fun(double *state, double dt, double *out_1385543446714745100) {
   out_1385543446714745100[0] = 1;
   out_1385543446714745100[1] = 0;
   out_1385543446714745100[2] = 0;
   out_1385543446714745100[3] = 0;
   out_1385543446714745100[4] = 0;
   out_1385543446714745100[5] = 0;
   out_1385543446714745100[6] = 0;
   out_1385543446714745100[7] = 0;
   out_1385543446714745100[8] = 0;
   out_1385543446714745100[9] = 1;
   out_1385543446714745100[10] = 0;
   out_1385543446714745100[11] = 0;
   out_1385543446714745100[12] = 0;
   out_1385543446714745100[13] = 0;
   out_1385543446714745100[14] = 0;
   out_1385543446714745100[15] = 0;
   out_1385543446714745100[16] = 0;
   out_1385543446714745100[17] = 0;
   out_1385543446714745100[18] = 1;
   out_1385543446714745100[19] = 0;
   out_1385543446714745100[20] = 0;
   out_1385543446714745100[21] = 0;
   out_1385543446714745100[22] = 0;
   out_1385543446714745100[23] = 0;
   out_1385543446714745100[24] = 0;
   out_1385543446714745100[25] = 0;
   out_1385543446714745100[26] = 0;
   out_1385543446714745100[27] = 1;
   out_1385543446714745100[28] = 0;
   out_1385543446714745100[29] = 0;
   out_1385543446714745100[30] = 0;
   out_1385543446714745100[31] = 0;
   out_1385543446714745100[32] = 0;
   out_1385543446714745100[33] = 0;
   out_1385543446714745100[34] = 0;
   out_1385543446714745100[35] = 0;
   out_1385543446714745100[36] = 1;
   out_1385543446714745100[37] = 0;
   out_1385543446714745100[38] = 0;
   out_1385543446714745100[39] = 0;
   out_1385543446714745100[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1385543446714745100[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1385543446714745100[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1385543446714745100[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1385543446714745100[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1385543446714745100[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1385543446714745100[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1385543446714745100[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1385543446714745100[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1385543446714745100[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1385543446714745100[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1385543446714745100[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1385543446714745100[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1385543446714745100[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1385543446714745100[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1385543446714745100[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1385543446714745100[56] = 0;
   out_1385543446714745100[57] = 0;
   out_1385543446714745100[58] = 0;
   out_1385543446714745100[59] = 0;
   out_1385543446714745100[60] = 0;
   out_1385543446714745100[61] = 0;
   out_1385543446714745100[62] = 0;
   out_1385543446714745100[63] = 1;
}
void h_25(double *state, double *unused, double *out_7483579869378611138) {
   out_7483579869378611138[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6605203896080581755) {
   out_6605203896080581755[0] = 0;
   out_6605203896080581755[1] = 0;
   out_6605203896080581755[2] = 0;
   out_6605203896080581755[3] = 0;
   out_6605203896080581755[4] = 0;
   out_6605203896080581755[5] = 0;
   out_6605203896080581755[6] = 1;
   out_6605203896080581755[7] = 0;
}
void h_24(double *state, double *unused, double *out_5260765027156920206) {
   out_5260765027156920206[0] = state[4];
   out_5260765027156920206[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7867592960276321887) {
   out_7867592960276321887[0] = 0;
   out_7867592960276321887[1] = 0;
   out_7867592960276321887[2] = 0;
   out_7867592960276321887[3] = 0;
   out_7867592960276321887[4] = 1;
   out_7867592960276321887[5] = 0;
   out_7867592960276321887[6] = 0;
   out_7867592960276321887[7] = 0;
   out_7867592960276321887[8] = 0;
   out_7867592960276321887[9] = 0;
   out_7867592960276321887[10] = 0;
   out_7867592960276321887[11] = 0;
   out_7867592960276321887[12] = 0;
   out_7867592960276321887[13] = 1;
   out_7867592960276321887[14] = 0;
   out_7867592960276321887[15] = 0;
}
void h_30(double *state, double *unused, double *out_6657927763590768282) {
   out_6657927763590768282[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3008695014868601525) {
   out_3008695014868601525[0] = 0;
   out_3008695014868601525[1] = 0;
   out_3008695014868601525[2] = 0;
   out_3008695014868601525[3] = 0;
   out_3008695014868601525[4] = 1;
   out_3008695014868601525[5] = 0;
   out_3008695014868601525[6] = 0;
   out_3008695014868601525[7] = 0;
}
void h_26(double *state, double *unused, double *out_1656254495929922637) {
   out_1656254495929922637[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6136947022689760588) {
   out_6136947022689760588[0] = 0;
   out_6136947022689760588[1] = 0;
   out_6136947022689760588[2] = 0;
   out_6136947022689760588[3] = 0;
   out_6136947022689760588[4] = 0;
   out_6136947022689760588[5] = 0;
   out_6136947022689760588[6] = 0;
   out_6136947022689760588[7] = 1;
}
void h_27(double *state, double *unused, double *out_2166854740064943338) {
   out_2166854740064943338[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4296277002705226837) {
   out_4296277002705226837[0] = 0;
   out_4296277002705226837[1] = 0;
   out_4296277002705226837[2] = 0;
   out_4296277002705226837[3] = 1;
   out_4296277002705226837[4] = 0;
   out_4296277002705226837[5] = 0;
   out_4296277002705226837[6] = 0;
   out_4296277002705226837[7] = 0;
}
void h_29(double *state, double *unused, double *out_2324041408416072483) {
   out_2324041408416072483[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4611892887621791475) {
   out_4611892887621791475[0] = 0;
   out_4611892887621791475[1] = 1;
   out_4611892887621791475[2] = 0;
   out_4611892887621791475[3] = 0;
   out_4611892887621791475[4] = 0;
   out_4611892887621791475[5] = 0;
   out_4611892887621791475[6] = 0;
   out_4611892887621791475[7] = 0;
}
void h_28(double *state, double *unused, double *out_3910161902706580420) {
   out_3910161902706580420[0] = state[5];
   out_3910161902706580420[1] = state[6];
}
void H_28(double *state, double *unused, double *out_7883068036422490090) {
   out_7883068036422490090[0] = 0;
   out_7883068036422490090[1] = 0;
   out_7883068036422490090[2] = 0;
   out_7883068036422490090[3] = 0;
   out_7883068036422490090[4] = 0;
   out_7883068036422490090[5] = 1;
   out_7883068036422490090[6] = 0;
   out_7883068036422490090[7] = 0;
   out_7883068036422490090[8] = 0;
   out_7883068036422490090[9] = 0;
   out_7883068036422490090[10] = 0;
   out_7883068036422490090[11] = 0;
   out_7883068036422490090[12] = 0;
   out_7883068036422490090[13] = 0;
   out_7883068036422490090[14] = 1;
   out_7883068036422490090[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
