#ifndef RAY_TRACE
#define RAY_TRACE

#include "camera.h"
#include "varr.h"


void trace(struct Camera *cam, double tMAX, int nhits, int ntracks,
            int (*target)(double, double *, double *, double *, void *),
            void *args);

void trace_single(int *status, int *iter, double *t1, double *xu1, 
                double t0, double *xu0, double tMAX, int n, struct varr *track,
                int (*target)(double, double *, double *, double *, void *),
                void *args);
void trace_xudot(double t, double *xu, void *args, double *xudot);

int target_eq_cart(double ta, double *xua, double *tb, double *xub, 
                    void *args);
int target_eq_sph(double ta, double *xua, double *tb, double *xub, void *args);

void trace_interpolateToCoordSurface(double *t, double *xu, 
                            double ta, double *xua, double tb, double *xub, 
                            int d, double C, void *args);

#endif
