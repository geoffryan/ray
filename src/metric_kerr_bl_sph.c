#include <math.h>
#include "par.h"
#include "metric.h"

// Kerr metric in Boyer-Lindquist coordinates.
// x = (t, r, theta, phi)

int metric_orientation_kerr_bl_sph()
{
    return SPH;
}

void  metric_g_kerr_bl_sph(double *g, double *x, void *args)
{
    int i;
    double r = x[1];
    double sint = sin(x[2]);
    double cost = cos(x[2]);
    double M = ((double *)args)[0];
    double a = ((double *)args)[1];

    double s2 = r*r + a*a*cost*cost;
    double del = r*r - 2*M*r + a*a;

    for(i=0; i<16; i++)
        g[i] = 0.0;
    
    g[0] = -1.0+2*M*r/s2;
    g[3] = -2*M*a*r*sint*sint/s2;
    g[5] = s2/del;
    g[10] = s2;
    g[12] = g[3];
    g[15] = (r*r + a*a*(1 - 2*M*r*sint*sint/s2)) * sint*sint;
}

void metric_dg_kerr_bl_sph(double *dg, double *x, void *args)
{
    int i;
    double r = x[1];
    double sint = sin(x[2]);
    double cost = cos(x[2]);
    double M = ((double *)args)[0];
    double a = ((double *)args)[1];

    double s2 = r*r + a*a*cost*cost;
    double del = r*r - 2*M*r + a*a;

    double ds2dr = 2*r;
    double ds2dth = -2*a*a*cost*sint;
    double ddeldr = 2*(r-M);

    for(i=0; i<64; i++)
        dg[i] = 0.0;
    
    dg[16*1+0] = 2*M*(s2 - r*ds2dr)/(s2*s2);
    dg[16*1+3] = -2*M*a*sint*sint * (s2 - r*ds2dr)/(s2*s2);
    dg[16*1+5] = (ds2dr*del - s2*ddeldr)/(del*del);
    dg[16*1+10] = ds2dr;
    dg[16*1+12] = dg[16*1+3];
    dg[16*1+15] = (2*r - 2*M*a*a*sint*sint*(s2-r*ds2dr)/(s2*s2)) * sint*sint;
    dg[16*2+0] = -2*M*r*ds2dth/(s2*s2);
    dg[16*2+3] = 2*M*a*r*(2*sint*cost*s2 - sint*sint*ds2dth)/(s2*s2);
    dg[16*2+5] = ds2dth/del;
    dg[16*2+10] = ds2dth;
    dg[16*2+12] = dg[16*2+3];
    dg[16*2+15] = 2*(r*r + a*a)*sint*cost
                - 2*M*r*a*a*sint*sint*sint*(4*cost*s2 - sint*ds2dth)/(s2*s2);
}

int metric_shadow_kerr_bl_sph(double *x, void *args)
{
    double M = ((double *)args)[0];
    double a = ((double *)args)[1];

    if(a > M)
        return 0;
    else if(x[1] < (1.0+1.0e-5)*M + sqrt((M-a)*(M+a)))
        return 1;
    return 0;
}

int metric_fix_domain_kerr_bl_sph(double *x, double *u, void *args)
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

