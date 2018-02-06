#include <math.h>
#include "par.h"
#include "metric.h"

// Schwarzschild metric in cartesian Kerr-Schild coordinates
// x = (t, x, y, z)
// args = (M)

int metric_orientation_schw_ks_cart()
{
    return CART;
}

void  metric_g_schw_ks_cart(double *g, double *X, void *args)
{
    double M = ((double *)args)[0];
    int i;
    for(i=0; i<16; i++)
        g[i] = 0.0;

    double x = X[1];
    double y = X[2];
    double z = X[3];

    double R = sqrt(x*x + y*y + z*z);
    double H = M/R;

    double k[4] = {1.0, x/R, y/R, z/R};

    int mu, nu;
    for(mu=0; mu<4; mu++)
        for(nu=0; nu<4; nu++)
        {
            if(mu == nu && mu == 0)
                g[4*mu+nu] = -1.0;
            else if(mu == nu)
                g[4*mu+nu] = 1.0;
            else
                g[4*mu+nu] = 0.0;
            g[4*mu+nu] += 2*H*k[mu]*k[nu];
        }
}

void metric_dg_schw_ks_cart(double *dg, double *X, void *args)
{
    int i, mu, nu;
    
    double x = X[1];
    double y = X[2];
    double z = X[3];

    double M = ((double *)args)[0];
    double R = sqrt(x*x + y*y + z*z);
    double H = 2*M/R;
    double k[4] = {1.0, x/R, y/R, z/R};
    double dkdXi[4] = {0.0, 0.0, 0.0, 0.0};
    double dXdXi[4] = {0.0, 0.0, 0.0, 0.0};

    for(mu=0; mu<16; mu++)
        dg[mu] = 0.0;

    for(i=1; i<4; i++)
    {
        for(mu=0; mu<4; mu++)
            dXdXi[mu] = 0.0;
        dXdXi[i] = 1.0;
        double dRdXi = X[i]/R;
        double dHdXi = -H*dRdXi/R;
        for(mu=1;mu<4; mu++)
            dkdXi[mu] = (R*dXdXi[mu] - X[mu]*dRdXi)/(R*R);

        for(mu=0; mu<4; mu++)
            for(nu=0; nu<=mu; nu++)
                dg[16*i+4*mu+nu] = dHdXi*k[mu]*k[nu]
                                    + H*(dkdXi[mu]*k[nu]+k[mu]*dkdXi[nu]);
        for(mu=0; mu<4; mu++)
            for(nu=mu+1; nu<4; nu++)
                dg[16*i+4*mu+nu] = dg[16*i+4*nu+mu];
    }
}

void metric_dg_schw_ks_nd_cart(double *dg, double *x, void *args)
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

int metric_shadow_schw_ks_cart(double *X, void *args)
{
    double M = ((double *)args)[0];
    double x = X[1];
    double y = X[2];
    double z = X[3];

    double R = sqrt(x*x + y*y + z*z);
    if(R < (2.0 + 1.0e-5)*M)
        return 1;
    return 0;
}

int metric_fix_domain_schw_ks_cart(double *x, double *u, void *args)
{
    return 0;
}

