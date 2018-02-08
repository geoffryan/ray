#ifndef RAY_METRIC
#define RAY_METRIC

#include <stdlib.h>
#include "par.h"

/*
 * metric_g(g, x, args):   Fills the instantiated array g with the values of
 *                          g_{\mu\nu} at position x^\mu and parameters 'args'
 *                          It is assumed 'g' has at least length 16, and x
 *                          has length at least 4.
 *
 * metric_dg(dg, x, args): Fills the instantiated array 'dg' with the values
 *                          of \partial_\mu g_{\nu\rho} at position x^\mu and 
 *                          parameters 'args'.  It is assumed 'dg' has length
 *                          at least 64, and x has length at least 4.  'dg' is 
 *                          filled with the derivative as the outer-most index.
 *                          ie. d_i g_[jk] == dg[16*i+4*j+k]
 */

enum{CART, SPH};

int (*metric_orientation)();
void (*metric_g)(double *, double *, void *);
void (*metric_dg)(double *, double *, void *);
int (*metric_shadow)(double *, void *);
int (*metric_fix_domain)(double *, double *, void *);
void (*metric_genPosFromSph)(double dist, double inc, double az, double *X);

int setup_metric(struct parList *pars);
void metric_genPosFromSph_cart(double dist, double inc, double az, double *X);
void metric_genPosFromSph_sph(double dist, double inc, double az, double *X);
void metric_ig(double *ig, double *x, double *g, void *args);
void metric_tetrad_euler(double *e, double *x, double *g, void *args);
void metric_tetrad(double *e, double *x, double *u, double *g, void *args);
void metric_null_ray(double *X, double tC, double pC, double *u, void *args);

int metric_orientation_flat_cart();
void metric_g_flat_cart(double *g, double *x, void *args);
void metric_dg_flat_cart(double *dg, double *x, void *args);
int metric_shadow_flat_cart(double *x, void *args);
int metric_fix_domain_flat_cart(double *x, double *u, void *args);

int metric_orientation_flat_sph();
void metric_g_flat_sph(double *g, double *x, void *args);
void metric_dg_flat_sph(double *dg, double *x, void *args);
int metric_shadow_flat_sph(double *x, void *args);
int metric_fix_domain_flat_sph(double *x, double *u, void *args);

int metric_orientation_schw_sc_sph();
void metric_g_schw_sc_sph(double *g, double *x, void *args);
void metric_dg_schw_sc_sph(double *dg, double *x, void *args);
int metric_shadow_schw_sc_sph(double *x, void *args);
int metric_fix_domain_schw_sc_sph(double *x, double *u, void *args);

int metric_orientation_schw_ks_sph();
void metric_g_schw_ks_sph(double *g, double *x, void *args);
void metric_dg_schw_ks_sph(double *dg, double *x, void *args);
int metric_shadow_schw_ks_sph(double *x, void *args);
int metric_fix_domain_schw_ks_sph(double *x, double *u, void *args);

int  metric_orientation_schw_ks_cart();
void metric_g_schw_ks_cart(double *g, double *x, void *args);
void metric_dg_schw_ks_cart(double *dg, double *x, void *args);
int metric_shadow_schw_ks_cart(double *x, void *args);
int metric_fix_domain_schw_ks_cart(double *x, double *u, void *args);

void metric_dg_schw_ks_nd_cart(double *dg, double *x, void *args);

int metric_orientation_kerr_bl_sph();
void metric_g_kerr_bl_sph(double *g, double *x, void *args);
void metric_dg_kerr_bl_sph(double *dg, double *x, void *args);
int metric_shadow_kerr_bl_sph(double *x, void *args);
int metric_fix_domain_kerr_bl_sph(double *x, double *u, void *args);

int  metric_orientation_bin_ks_cart();
void  metric_g_bin_ks_cart(double *g, double *x, void *args);
void metric_dg_bin_ks_cart(double *dg, double *x, void *args);
int metric_shadow_bin_ks_cart(double *x, void *args);
int metric_fix_domain_bin_ks_cart(double *x, double *u, void *args);

#endif

