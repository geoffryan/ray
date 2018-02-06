#include <math.h>
#include "par.h"
#include "metric.h"

// Schwarzschild metric in spherical Kerr-Schild coordinates.
// x = (t, r, theta, phi)

int metric_orientation_schw_ks_sph()
{
    return SPH;
}

void  metric_g_schw_ks_sph(double *g, double *x, void *args)
{
    int i;
    double r = x[1];
    double sint = sin(x[2]);
    double M = ((double *)args)[0];

    for(i=0; i<16; i++)
        g[i] = 0.0;
    
    g[0] = -1.0+2*M/r;
    g[1] = 2*M/r;
    g[4] = g[1];
    g[5] = 1.0 + 2*M/r;
    g[10] = r*r;
    g[15] = r*r*sint*sint;
}

void metric_dg_schw_ks_sph(double *dg, double *x, void *args)
{
    int i;
    for(i=0; i<64; i++)
        dg[i] = 0.0;

    double M = ((double *)args)[0];
    double r = x[1];
    double sint = sin(x[2]);
    double cost = cos(x[2]);

    dg[16*1 + 0] = -2*M/(r*r);
    dg[16*1 + 1] = -2*M/(r*r);
    dg[16*1 + 4] = -2*M/(r*r);
    dg[16*1 + 5] = -2*M/(r*r);
    dg[16*1 + 10] = 2.0*r;
    dg[16*1 + 15] = 2.0*r*sint*sint;
    dg[16*2 + 15] = 2.0*r*r*sint*cost;
}

int metric_shadow_schw_ks_sph(double *x, void *args)
{
    double M = ((double *)args)[0];

    if(x[1] < (2.0 + 1.0e-5)*M)
        return 1;
    return 0;
}

int metric_fix_domain_schw_ks_sph(double *x, double *u, void *args)
{
    if(x[2] < 0.0)
    {
        x[2] = -x[2];
        x[3] += M_PI;
        u[2] = -u[2];
        return 1;
    }

    return 0;
}

