#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "metric.h"
#include "ode.h"
#include "output.h"
#include "surface.h"
#include "trace.h"

void trace(struct Camera *cam, double tMAX, int maxhits, int ntracks,
            void *args, char *filename, int fancy_printing)
{
    int mu, nu, i, n;
    double g[16];

    int N = cam->N;
    srand(27181);
    rand();
    int tids[ntracks];
    for(n=0; n<ntracks; n++)
        tids[n] = rand() % N;

    maxhits = maxhits > 1 ? maxhits : 1;
    int buf_width = 2 + 9*(maxhits + 1);

    double *map = (double *) malloc(buf_width*N*sizeof(double));
    int *nhits = (int *) malloc(N*sizeof(int));
    struct varr track = VARR_DEFAULT;
    varr_init(&track, 1800);

    metric_g(g, cam->X, args);

    int Np = cam->Np;
    int Nt = cam->Nt;
    double aspect = fabs((cam->phiMax-cam->phiMin)
                            /(cam->thetaMax-cam->thetaMin));

    int Ny = 16;
    int Nx = 2*Ny*aspect*1.01;
    char line[Nx+1];
    double avgs[Nx];
    double nx = ((double)Np) / Nx;
    double ny = ((double)Nt) / Ny;

    for(i=0; i<Nx; i++)
        line[i] = ' ';
    line[Nx] = '\0';
    for(i=0; i<Nx; i++)
        avgs[i] = 0.0;

    if(fancy_printing)
        printf("\r%s %.1lf%%", line, 0.0);
    else
        printf("%.1lf%%", 0.0);

    for(n=0; n<N; n++)
    {
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

        /*
        printf("   x: %.6lf %.6lf %.6lf %.6lf\n",
                    xu0[0], xu0[1], xu0[2], xu0[3]);
        printf("   u: %.6lf %.6lf %.6lf %.6lf\n",
                    xu0[4], xu0[5], xu0[6], xu0[7]);
        */

        double t0, tHits[maxhits], xuHits[8*maxhits];
        int status, iter, id;
        id = -1;
        for(i=0; i<ntracks; i++)
            if(tids[i] == n)
                id = n;
        t0 = 0.0;

        status = trace_single(&iter, &(nhits[n]), tHits, xuHits, t0, xu0, tMAX,
                                maxhits, id, &track, args);


        //printf("    %d iterations - t: %.6lf\n", iter, tHits[nhits[n]-1]);

        map[buf_width * n + 0] = cam->thC[2*n];
        map[buf_width * n + 1] = cam->thC[2*n+1];
        map[buf_width * n + 2] = t0;
        for(i=0; i<8; i++)
            map[buf_width*n+i+3] = xu0[i];
        int j;
        for(j=0; j<maxhits; j++)
        {
            map[buf_width * n + 9*j + 11] = tHits[j];
            for(i=0; i<8; i++)
                map[buf_width*n + 9*j + i + 12] = xuHits[8*j+i];
        }

        if(id >= 0)
            output_track_h5(&track, id, filename);
        
        varr_clear(&track);

        avgs[(int)((n%Np) / nx)] += ((double)status) / (nx*ny);
        
        if(fancy_printing
            && (((int) ((float)1000*n)/N) != ((int) ((float)1000*(n+1))/N)
                || (n+1)%Np == 0 ||n == N-1))
        {
            for(i=0; i<Nx; i++)
            {
                if(avgs[i] < 0.0)
                    line[i] = 'o';
                else if(avgs[i] < 1.0)
                    line[i] = '*';
                else
                    line[i] = '+';
            }
            if(((int) (((n+1)/Np)/ny)) != ((int) ((n/Np)/ny)) && n<N-1)
            {
                printf("\r%s       \n", line);
                for(i=0; i<Nx; i++)
                {
                    line[i] = ' ';
                    avgs[i] = 0.0;
                }
            }
            printf("\r%s %.1lf%%", line, (100.0*(n+1))/N);
            fflush(stdout);
        }
        else if( ((int) ((float)1000*n)/N) != ((int) ((float)1000*(n+1))/N) )
        {
            printf("\r%.1lf%%", (100.0*(n+1))/N);
            fflush(stdout);
        }
    }
    printf("\n");

    output_map_h5(map, nhits, N, buf_width, filename);

    free(map);
    printf("Track varr had size %d (%d steps)\n", track.size, track.size/9);
    varr_free(&track);
}

int trace_single(int *iter, int *nhits, double *tHits, double *xuHits, 
                double t0, double *xu0, double tMAX, int maxhits, 
                int id, struct varr *track,
                void *args)
{
        double xu[8];
        double xu1[8];
        int mu;

        for(mu=0; mu<8; mu++)
        {
            xu[mu] = xu0[mu];
            xu1[mu] = xu0[mu];
        }

        double t = t0;
        double dt = -1.0e-6;
        double dt0;

        int hit_status = 0;

        if(id >= 0)
        {
            varr_clear(track);
            varr_append(track, t);
            varr_append_chunk(track, xu0, 8);
        }

        int i = 0;
        int ihit = 0;
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
            
            hit_status = trace_target(t-dt0, xu, t, xu1, &(tHits[ihit]),
                                        &(xuHits[8*ihit]), args);
            if(hit_status != 0)
                ihit++;

            if(id >= 0)
            {
                if(hit_status != 0)
                {
                    varr_append(track, tHits[ihit-1]);
                    varr_append_chunk(track, &(xuHits[8*(ihit-1)]), 8);
                }
                else
                {
                    varr_append(track, t);
                    varr_append_chunk(track, xu1, 8);
                }
            }

            //printf("        target - status: %.6d\n", status);

            if(hit_status < 0 || ihit == maxhits)
                break;
        }

        if(fabs(t) > fabs(tMAX) && ihit < maxhits)
        {
            tHits[ihit] = t;
            for(mu=0; mu<8; mu++)
                xuHits[8*ihit+mu] = xu1[mu];
            ihit++;
        }

        *iter = i;
        *nhits = ihit;

        return hit_status;
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


int trace_target(double ta, double *xua, double tb, double *xub, 
                    double *t, double *xu, void *args)
{
    int status = 0;

    if(metric_shadow(xub, args))
    {
        status = -1;
        *t = tb;
        int mu;
        for(mu=0; mu<8; mu++)
            xu[mu] = xub[mu];
    }
    else
    {
        double fa = surface_f(ta, xua, args);
        double fb = surface_f(tb, xub, args);

        if(fb*fa < 0)
        {
            status = 1;
            
            trace_interpolateToSurface(t, xu, ta, xua, tb, xub, args);
        }
    }

    return status;
}

void trace_interpolateToSurface(double *t, double *xu, double ta, double *xua,
                                double tb, double *xub, void *args)
{
    double atol = 1.0e-8;

    double xudota[8], xudotb[8];
    trace_xudot(ta, xua, args, xudota);
    trace_xudot(tb, xub, args, xudotb);

    // Calculate interpolating polynomial for target variable
    double fa = surface_f(ta, xua, args);
    double fb = surface_f(tb, xub, args);
    double ma = surface_df(ta, xua, xudota, args);
    double mb = surface_df(tb, xub, xudotb, args);

    double c0, c1, c2, c3;
    c0 = fa;
    c1 = (tb-ta)*ma;
    c2 = -3*fa - 2*(tb-ta)*ma + 3*fb - (tb-ta)*mb;
    c3 = 2*fa + (tb-ta)*ma - 2*fb + (tb-ta)*mb;
  
    /*
    printf("ta: %.14lg\n", ta);
    printf("xua: %.14lg %.14lg %.14lg %.14lg\n",
            xua[0], xua[1], xua[2], xua[3]);
    printf("xudota: %.14lg %.14lg %.14lg %.14lg\n",
            xudota[0], xudota[1], xudota[2], xudota[3]);
    printf("tb: %.14lg\n", tb);
    printf("xub: %.14lg %.14lg %.14lg %.14lg\n",
            xub[0], xub[1], xub[2], xub[3]);
    printf("xudotb: %.14lg %.14lg %.14lg %.14lg\n",
            xudotb[0], xudotb[1], xudotb[2], xudotb[3]);
    */

    //Find intersection of interpolant and surface
    double s, ds, f, df;
    s = 0.5;
    do
    {
        f = c0 + s*(c1 + s*(c2 + s*c3));
        df = c1 + 2*s*c2 + 3*s*s*c3;
        ds = -f/df;
        //printf("   %.3lf %.3lf %.3lf\n", s, f, df);
        s += ds;
    }
    while(fabs(ds) > atol && fabs(f) > atol);

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
    printf("s: %.14lg\n", s);
    printf("t: %.14lg\n", *t);
    printf("xu: %.14lg %.14lg %.14lg %.14lg\n",
            xu[0], xu[1], xu[2], xu[3]);
    */
}
