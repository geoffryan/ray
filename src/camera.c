#include <stdlib.h>
#include <stdio.h>
#include "metric.h"
#include "camera.h"

int setup_camera(struct Camera *cam, struct parList *par)
{
    int meshType = 0;
    int Nt = par->N1;
    int Np = par->N2;
    int N = Nt*Np;

    cam->meshType = meshType;
    cam->N = N;
    cam->Nt = Nt;
    cam->Np = Np;

    double ta = par->X1a * M_PI;
    double tb = par->X1b * M_PI;
    double pa = par->X2a * M_PI;
    double pb = par->X2b * M_PI;

    cam->thetaMax = tb;
    cam->thetaMin = ta;
    cam->phiMax = pb;
    cam->phiMin = pa;

    double *t = (double *)malloc(Nt * sizeof(double));
    double *p = (double *)malloc(Np * sizeof(double));

    int i;
    for(i=0; i<Nt; i++)
        t[i] = ta + (i+0.5)*(tb-ta)/Nt;
    for(i=0; i<Np; i++)
        p[i] = pa + (i+0.5)*(pb-pa)/Np;

    cam->thC = (double *) malloc(2 * N * sizeof(double));

    for(i=0; i<N; i++)
    {
        cam->thC[2*i] = t[i/Np];
        cam->thC[2*i+1] = p[i%Np];
    }

    free(t);
    free(p);

    return 0;
}

void camera_print(struct Camera *cam)
{
    printf("X: (%.6e, %.6e, %.6e, %.6e)\n", cam->X[0], cam->X[1], cam->X[2],
                                            cam->X[3]);
    printf("U: (%.6e, %.6e, %.6e, %.6e)\n", cam->U[0], cam->U[1], cam->U[2],
                                            cam->U[3]);
    printf("e: (%.6e, %.6e, %.6e, %.6e)\n", cam->e[0], cam->e[1], cam->e[2],
                                            cam->e[3]);
    printf("   (%.6e, %.6e, %.6e, %.6e)\n", cam->e[4], cam->e[5], cam->e[6],
                                            cam->e[7]);
    printf("   (%.6e, %.6e, %.6e, %.6e)\n", cam->e[8], cam->e[9], cam->e[10],
                                            cam->e[11]);
    printf("   (%.6e, %.6e, %.6e, %.6e)\n", cam->e[12], cam->e[13], cam->e[14],
                                            cam->e[15]);
    printf("N: %d\n", cam->N);
}

void camera_set_position(struct Camera *cam, double *x, double *u, void *args)
{
    double g[16];
    double e[16];

    metric_g(g, x, args);
    metric_tetrad(e, x, u, g, args);

    int mu;
    for(mu=0; mu<4; mu++)
    {
        cam->X[mu] = x[mu];
        cam->U[mu] = u[mu];
    }
    for(mu=0; mu<16; mu++)
        cam->e[mu] = e[mu];
}

void camera_get_ray(struct Camera *cam, double *k, int id)
{
    double t = cam->thC[2*id];
    double p = cam->thC[2*id+1];

    double cost = cos(t);
    double sint = sin(t);
    double cosp = cos(p);
    double sinp = sin(p);

    double kc[4];
    // In general, the desired orientation is for a right-handed system, with
    // "forward" (phi=0,theta=pi/2) towards the origin and "down" (theta=pi)
    // towards -z.
    if(metric_orientation() == SPH)
    {
        kc[0] = 1.0;
        kc[1] = sint*cosp;
        kc[2] = cost;
        kc[3] = sint*sinp;
    }
    else
    {
        double x = cam->X[1];
        double y = cam->X[2];
        double z = cam->X[3];
        double r = sqrt(x*x+y*y);
        double R = sqrt(x*x+y*y+z*z);
        double cosi = z/R;
        double sini = r/R;
        double cosa = x/r;
        double sina = y/r;

        // "Local" frame
        double kc0[4] = {1.0, -sint*cosp, -sint*sinp, -cost};
        //Pitch up to align global z
        double kc1[4] = {kc0[0], sini*kc0[1]+cosi*kc0[3],
                            kc0[2], -cosi*kc0[1]+sini*kc0[3]};
        //Rotate (yaw) to align x,y
        double kc2[4] = {kc1[0], -cosa*kc1[1]+sina*kc1[2],
                            -sina*kc1[1]-cosa*kc1[2], kc1[3]};

        kc[0] = kc2[0];
        kc[1] = kc2[1];
        kc[2] = kc2[2];
        kc[3] = kc2[3];
    }

    int mu, i;
    for(mu=0; mu<4; mu++)
    {
        k[mu] = 0.0;
        for(i=0; i<4; i++)
            k[mu] += cam->e[4*mu+i]*kc[i];
    }
}

void free_camera(struct Camera *cam)
{
    free(cam->thC);
}
