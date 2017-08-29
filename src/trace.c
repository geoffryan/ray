#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "metric.h"
#include "ode.h"
#include "trace.h"

void trace(struct Camera *cam, double tMAX, int track_thin,
            int target(double, double *, double *, double *, void *),
            void *args)
{
    int mu, nu, i, n;
    double g[16];
    FILE *f;

    int N = cam->N;

    double *map = (double *) malloc(19*N*sizeof(double));

    metric_g(g, cam->X, args);

    for(n=0; n<N; n++)
    {
        printf("Integrating ray %d\n", n);
        double xu0[8];
        for(mu=0; mu<4; mu++)
            xu0[mu] = cam->X[mu];
        double k[4];
        camera_get_ray(cam, k, n);
        for(mu=0; mu<4; mu++)
        {
            xu0[4+mu] = 0.0;
            for(nu=0; nu<4; nu++)
                xu0[4+mu] += g[4*mu+nu]*k[nu];
        }

        printf("   x: %.6lf %.6lf %.6lf %.6lf\n",
                    xu0[0], xu0[1], xu0[2], xu0[3]);
        printf("   u: %.6lf %.6lf %.6lf %.6lf\n",
                    xu0[4], xu0[5], xu0[6], xu0[7]);

        double xu[8];
        double xu1[8];
        for(mu=0; mu<8; mu++)
        {
            xu[mu] = xu0[mu];
            xu1[mu] = xu0[mu];
        }

        double t = 0;
        double dt = -1.0e-6;
        double dt0;
        int status;

        if(n % track_thin == 0)
        {
            char fname[256];
            sprintf(fname, "track_%04d.txt", n);

            f = fopen(fname, "w");
            fprintf(f, "%.12lg", t);
            for(mu=0; mu<8; mu++)
                fprintf(f, " %.12lg", xu0[mu]);
            fprintf(f, "\n");
        }

        i=0;
        while(fabs(t) < fabs(tMAX))
        {
            //printf("      %d: %.6lf\n", i, t);
            for(mu=0; mu<8; mu++)
                xu[mu] = xu1[mu];

            dt0 = dopr54(t, xu1, &dt, 8, args, trace_xudot);
            metric_fix_domain(&(xu[0]), &(xu[4]), args);
            //printf("        step - dt: %.6lf dt0: %.6lf\n", dt, dt0);
            t += dt0;
            i++;
            
            if(n % track_thin == 0)
            {
                fprintf(f, "%.12lg", t);
                for(mu=0; mu<8; mu++)
                    fprintf(f, " %.12lg", xu1[mu]);
                fprintf(f, "\n");
            }

            status = target(t-dt0, xu, &t, xu1, args);
            //printf("        target - status: %.6d\n", status);

            if(status != 0)
                break;

        }

        if(n % track_thin == 0)
        {
            fclose(f);
        }

        printf("    %d iterations - t: %.6lf\n", i, t);

        map[19*n + 0] = cam->thC[2*n];
        map[19*n + 1] = cam->thC[2*n+1];
        map[19*n + 2] = t;
        for(i=0; i<8; i++)
            map[19*n+i+3] = xu0[i];
        for(i=0; i<8; i++)
            map[19*n+i+11] = xu1[i];
    }

    char fname[] = "map.txt";
    f = fopen(fname, "w");
    for(n=0; n<N; n++)
    {
        fprintf(f, "%.6lf", map[19*n]);
        for(i=1; i<19; i++)
            fprintf(f, " %.6lf", map[19*n+i]);
        fprintf(f, "\n");
    }
    fclose(f);

    free(map);
}

void trace_xudot(double t, double *xu, void *args, double *xudot)
{
    double *x = &(xu[0]);
    double *u = &(xu[4]);
    double *xdot = &(xudot[0]);
    double *udot = &(xudot[4]);

    int mu, nu, si;
    double g[16], ig[16], dg[64], uu[4];
    metric_g(g, x, args);
    metric_ig(ig, x, g, args);
    metric_dg(dg, x, args);

    for(mu=0; mu<4; mu++)
    {
        uu[mu] = 0.0;
        for(nu=0; nu<4; nu++)
            uu[mu] += ig[4*mu+nu]*u[nu];
    }

    for(mu=0; mu<4; mu++)
        xdot[mu] = uu[mu];

    for(mu=0; mu<4; mu++)
    {
        udot[mu] = 0;
        for(nu=0; nu<4; nu++)
            for(si=0; si<4; si++)
                udot[mu] += 0.5*uu[nu]*uu[si]*dg[16*mu+4*nu+si];
    }
}

int target_eq_cart(double ta, double *xua, double *tb, double *xub, void *args)
{
    int status = 0;

    if(metric_shadow(xub, args))
        status = -1;
    else
    {
        double za = xua[3];
        double zb = xub[3];

        if(zb*za < 0)
            status = 1;
    }

    return status;
}

int target_eq_sph(double ta, double *xua, double *tb, double *xub, void *args)
{
    int status = 0;

    if(metric_shadow(xub, args))
        status = -1;
    else
    {
        double tha = xua[2];
        double thb = xub[2];

        if((tha-0.5*M_PI)*(thb-0.5*M_PI) < 0)
        {
            status = 1;
            double t, xu[8];
            trace_interpolateToCoordSurface(&t, xu, ta, xua, *tb, xub,
                                    2, 0.5*M_PI, args);
            *tb = t;
            int mu;
            for(mu=0; mu<8; mu++)
                xub[mu] = xu[mu];
        }
    }

    return status;
}

void trace_interpolateToCoordSurface(double *t, double *xu, 
                            double ta, double *xua, double tb, double *xub, 
                            int d, double C, void *args)
{
    double atol = 1.0e-8;
    double rtol = 1.0e-8;

    double xudota[8], xudotb[8];
    trace_xudot(ta, xua, args, xudota);
    trace_xudot(tb, xub, args, xudotb);

    // Calculate interpolating polynomial for target variable
    double fa = xua[d];
    double fb = xub[d];
    double ma = xudota[d];
    double mb = xudotb[d];

    double c0, c1, c2, c3;
    c0 = fa;
    c1 = (tb-ta)*ma;
    c2 = -3*fa - 2*(tb-ta)*ma + 3*fb - (tb-ta)*mb;
    c3 = 2*fa + (tb-ta)*ma - 2*fb + (tb-ta)*mb;
  
    /*
    printf("ta: %.14lg\n", ta);
    printf("xua: %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg\n",
            xua[0], xua[1], xua[2], xua[3], xua[4], xua[5], xua[6], xua[7]);
    printf("xudota: %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg\n",
            xudota[0], xudota[1], xudota[2], xudota[3], xudota[4], xudota[5],
            xudota[6], xudota[7]);
    printf("tb: %.14lg\n", tb);
    printf("xub: %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg\n",
            xub[0], xub[1], xub[2], xub[3], xub[4], xub[5], xub[6], xub[7]);
    printf("xudotb: %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg %.14lg\n",
            xudotb[0], xudotb[1], xudotb[2], xudotb[3], xudotb[4], xudotb[5],
            xudotb[6], xudotb[7]);
    */

    //Find intersection of interpolant and surface
    double s, ds, f, df;
    s = 0.5;
    do
    {
        f = c0 + s*(c1 + s*(c2 + s*c3)) - C;
        df = c1 + 2*s*c2 + 3*s*s*c3;
        ds = -f/df;
        //printf("   %.3lf %.3lf %.3lf\n", s, f, df);
        s += ds;
    }
    while(fabs(ds) > atol && fabs(f-C) > atol + fabs(C)*rtol);

    //Calculate all variables at surface location
    double h00 = (1+2*s)*(1-s)*(1-s);
    double h10 = s*(1-s)*(1-s);
    double h01 = s*s*(3-2*s);
    double h11 = s*s*(s-1);

    int mu;
    for(mu=0; mu<8; mu++)
    {
        xu[mu] = h00*xua[mu] + h10*(tb-ta)*xudota[mu] + h01*xub[mu]
                    + h11*(tb-ta)*xudotb[mu];
    }
    *t = (1-s)*ta + s*tb;

    /*
    printf("s: %.6lf\n", s);
    printf("xu: %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",
            xu[0], xu[1], xu[2], xu[3], xu[4], xu[5], xu[6], xu[7]);
    */
}
