#include <math.h>
#include "par.h"
#include "metric.h"

// Flat space metric in spherical polar coordinates.
// x = (t, r, theta, phi)

void  metric_g_flat_sph(double *g, double *x, void *args)
{
    int i;
    double sint = sin(x[2]);
    for(i=0; i<16; i++)
        g[i] = 0.0;
    g[0] = -1.0;
    g[5] = 1.0;
    g[10] = x[1]*x[1];
    g[15] = x[1]*x[1]*sint*sint;
}

void metric_dg_flat_sph(double *dg, double *x, void *args)
{
    int i;
    for(i=0; i<64; i++)
        dg[i] = 0.0;
    
    double sint = sin(x[2]);
    double cost = cos(x[2]);

    dg[16*1 + 10] = 2.0*x[1];
    dg[16*1 + 15] = 2.0*x[1]*sint*sint;
    dg[16*2 + 15] = 2.0*cost*x[1]*x[1]*sint;
}

int metric_shadow_flat_sph(double *x, void *args)
{
    return 0;
}

int metric_fix_domain_flat_sph(double *x, double *u, void *args)
{
    if(x[2] < 0.0)
    {
        x[2] = -x[2];
        x[3] += M_PI;
        u[2] = -u[2];
        return 1;
    }
    if(x[2] > M_PI)
    {
        x[2] = 2*M_PI - x[2];
        x[3] += M_PI;
        u[2] = -u[2];
    }

    return 0;
}

