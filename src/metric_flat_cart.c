#include <math.h>
#include "par.h"
#include "metric.h"

// Flat space metric in cartesian coordinates.
// x = (t, x, y, z)

int metric_orientation_flat_cart()
{
    return CART;
}

void  metric_g_flat_cart(double *g, double *x, void *args)
{
    int i;
    for(i=0; i<16; i++)
        g[i] = 0.0;

    g[0] = -1.0;
    g[5] = 1.0;
    g[10] = 1.0;
    g[15] = 1.0;
}

void  metric_ig_flat_cart(double *g, double *x, void *args)
{
    int i;
    for(i=0; i<16; i++)
        g[i] = 0.0;

    g[0] = -1.0;
    g[5] = 1.0;
    g[10] = 1.0;
    g[15] = 1.0;
}

void metric_dg_flat_cart(double *dg, double *x, void *args)
{
    int i;
    for(i=0; i<64; i++)
        dg[i] = 0.0;
}

int metric_shadow_flat_cart(double *x, void *args)
{
    return 0;
}

int metric_fix_domain_flat_cart(double *x, double *u, void *args)
{
    return 0;
}

