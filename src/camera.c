#include <stdlib.h>
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

    double ta = par->X1a * M_PI;
    double tb = par->X1b * M_PI;
    double pa = par->X2a * M_PI;
    double pb = par->X2b * M_PI;

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

    // For Cartesian setups
    //double kc[4] = {1.0, -sint*cosp, -sint*sinp, -cost};

    //For Spherical setups ("down" is towards equator, "forwards" is toward 
    //origin
    double kc[4] = {1.0, sint*cosp, cost, sint*sinp};

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
