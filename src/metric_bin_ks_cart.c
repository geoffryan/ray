#include <math.h>
#include "par.h"
#include "metric.h"

// Approximate binary BH metric in cartesian coordinates.
// x = (t, x, y, z)
// args = (M, q, a, om)

int metric_orientation_bin_ks_cart()
{
    return CART;
}

void calc_pars(double *M1, double *M2, double *a1, double *a2, double *om,
                void *args)
{
    double M = ((double *)args)[0];
    double q = ((double *)args)[1];
    double a = ((double *)args)[2];
    *om = ((double *)args)[3];

    *M1 = M/(1+q);
    *M2 = q*M/(1+q);
    *a1 = q*a/(1+q);
    *a2 = a/(1+q);
}

void  metric_g_bin_ks_cart(double *g, double *X, void *args)
{
    int i,j;
    for(i=0; i<16; i++)
        g[i] = 0.0;

    double M1, M2, a1, a2, om;
    calc_pars(&M1, &M2, &a1, &a2, &om, args);

    double t = X[0];
    double x = X[1];
    double y = X[2];
    double z = X[3];

    double phiBH = om*t;
    double x1 = a1*cos(phiBH);
    double y1 = a1*sin(phiBH);
    double x2 = -a2*cos(phiBH);
    double y2 = -a2*sin(phiBH);

    double R1 = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z);
    double R2 = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z);
    double H1 = M1/R1;
    double H2 = M2/R2;

    double k1[3] = {(x1-x)/R1, (y1-y)/R1, -z/R1};
    double k2[3] = {(x2-x)/R2, (y2-y)/R2, -z/R2};

    double al2 = 1.0 / ((1+2*H1)*(1+2*H2));
    double be[3];
    be[0] = -2*H1*k1[0]/(1+2*H1) - 2*H2*k2[0]/(1+2*H2);
    be[1] = -2*H1*k1[1]/(1+2*H1) - 2*H2*k2[1]/(1+2*H2);
    be[2] = -2*H1*k1[2]/(1+2*H1) - 2*H2*k2[2]/(1+2*H2);

    double gam[3][3];
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
        {
            if(i == j)
                gam[i][j] = 1.0;
            else
                gam[i][j] = 0.0;
            gam[i][j] += 2*H1*k1[i]*k1[j] + 2*H2*k2[i]*k2[j];
        }

    double b2, bed[3];
    b2 = 0.0;
    for(i=0; i<3; i++)
    {
        bed[i] = 0.0;
        for(j=0; j<3; j++)
            bed[i] += gam[i][j]*be[j];
        b2 += be[i]*bed[i];
    }

    g[0] = -al2 + b2;
    for(i=0; i<3; i++)
    {
        g[4*0+i+1] = bed[i];
        g[4*(i+1)+0] = bed[i];
        for(j=0; j<3; j++)
            g[4*(i+1)+j+1] = gam[i][j];
    }
}

void metric_dg_bin_ks_cart(double *dg, double *x, void *args)
{
    double dx = 1.0e-6;
    int mu;
    for(mu=0; mu<4; mu++)
    {
        double xp[4] = {x[0], x[1], x[2], x[3]};
        double xm[4] = {x[0], x[1], x[2], x[3]};
        xp[mu] += dx;
        xm[mu] -= dx;
        double gp[16], gm[16];
        metric_g(gp, xp, args);
        metric_g(gm, xm, args);

        int nu;
        double i2dx = 1.0/(2*dx);
        for(nu=0; nu<16; nu++)
            dg[16*mu + nu] = i2dx*(gp[nu]-gm[nu]);
    }
}

int metric_shadow_bin_ks_cart(double *X, void *args)
{
    double M1, M2, a1, a2, om;
    calc_pars(&M1, &M2, &a1, &a2, &om, args);

    double t = X[0];
    double x = X[1];
    double y = X[2];
    double z = X[3];

    double phiBH = om*t;
    double x1 = a1*cos(phiBH);
    double y1 = a1*sin(phiBH);
    double x2 = -a2*cos(phiBH);
    double y2 = -a2*sin(phiBH);

    double R1 = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z);
    double R2 = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z);

    if(R1 < (2.0 + 1.0e-5)*M1 || R2 < (2.0 + 1.0e-5)*M2)
        return 1;
    return 0;
}

int metric_fix_domain_bin_ks_cart(double *x, double *u, void *args)
{
    return 0;
}

