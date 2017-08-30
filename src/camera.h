#ifndef RAY_CAMERA
#define RAY_CAMERA

#include <math.h>
#include <stdlib.h>
#include "par.h"

struct Camera
{
    double X[4];
    double U[4];
    double e[16];

    int meshType;
    int N;

    double *thC;
};

const static struct Camera CAMERA_DEFAULT = {
    .X = {0., 0., 0., 0.},
    .U = {1., 0., 0., 0.},
    .e = {1.,0.,0.,0., 0.,1.,0.,0., 0.,0.,1.,0., 0.,0.,0.,1.},
    .meshType = 0,
    .N = 0,
    .thC = NULL,
};

int setup_camera(struct Camera *cam, struct parList *par);
int camera_initialize(struct Camera *cam, struct parList *par);
void camera_set_position(struct Camera *cam, double *x, double *u, void *args);
void camera_get_ray(struct Camera *cam, double *k, int n);
void free_camera(struct Camera *cam);

#endif
