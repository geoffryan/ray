#ifndef RAY_TRACE
#define RAY_TRACE

#include "camera.h"
#include "varr.h"

void trace(struct Camera *cam, double tMAX, int nhits, int ntracks, 
            void *args, char *filename);

int trace_single(int *status, int *iter, double *t1, double *xu1, 
                double t0, double *xu0, double tMAX, int ntracks, 
                int n, struct varr *track, void *args);
void trace_xudot(double t, double *xu, void *args, double *xudot);

int trace_target(double ta, double *xua, double tb, double *xub,
                    double *t, double *xu, void *args);
void trace_interpolateToSurface(double *t, double *xu, double ta, double *xua,
                                double tb, double *xub, void *args);

#endif
