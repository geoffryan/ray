#include <math.h>
#include <stdio.h>
#include "metric.h"
#include "surface.h"

int setup_surface(struct parList *pars)
{
    int choice = pars->surface;
    int err = 0;

    if(choice == 0)
    {
        surface_f = &surface_equator_f;
        surface_df = &surface_equator_df;
    }
    else if(choice == 1)
    {
        surface_f = &surface_sky_f;
        surface_df = &surface_sky_df;
    }
    else
    {
        printf("Bad Surface choice: %d\n", choice);
        err = 1;
    }

    return err;
}

double surface_equator_f(double t, double *xu, void *args)
{
    if(metric_orientation() == SPH)
        return 0.5*M_PI - xu[2];
    else
        return xu[3];
}

double surface_equator_df(double t, double *xu, double *xudot, void *args)
{
    if(metric_orientation() == SPH)
        return -xudot[2];
    else
        return xudot[3];
}

double surface_sky_f(double t, double *xu, void *args)
{
    double Rsky = 1.0e5;

    if(metric_orientation() == SPH)
        return xu[1] - Rsky;
    else
        return xu[1]*xu[1]+xu[2]*xu[2]+xu[3]*xu[3] - Rsky*Rsky;

}

double surface_sky_df(double t, double *xu, double *xudot, void *args)
{
    if(metric_orientation() == SPH)
        return xudot[1];
    else
        return 2*(xu[1]*xudot[1]+xu[2]*xudot[2]+xu[3]*xudot[3]);
}

