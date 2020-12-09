/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5091353622417362108);
void inv_err_fun(double *nom_x, double *true_x, double *out_3347397060309790859);
void H_mod_fun(double *state, double *out_8038179413808804489);
void f_fun(double *state, double dt, double *out_6958450483037231644);
void F_fun(double *state, double dt, double *out_1385543446714745100);
void h_25(double *state, double *unused, double *out_7483579869378611138);
void H_25(double *state, double *unused, double *out_6605203896080581755);
void h_24(double *state, double *unused, double *out_5260765027156920206);
void H_24(double *state, double *unused, double *out_7867592960276321887);
void h_30(double *state, double *unused, double *out_6657927763590768282);
void H_30(double *state, double *unused, double *out_3008695014868601525);
void h_26(double *state, double *unused, double *out_1656254495929922637);
void H_26(double *state, double *unused, double *out_6136947022689760588);
void h_27(double *state, double *unused, double *out_2166854740064943338);
void H_27(double *state, double *unused, double *out_4296277002705226837);
void h_29(double *state, double *unused, double *out_2324041408416072483);
void H_29(double *state, double *unused, double *out_4611892887621791475);
void h_28(double *state, double *unused, double *out_3910161902706580420);
void H_28(double *state, double *unused, double *out_7883068036422490090);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
