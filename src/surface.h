#ifndef RAY_SURFACE
#define RAY_SURFACE

double (*surface_f)(double, double *, void *);
double (*surface_df)(double, double *, double *, void *);

int setup_surface(struct parList *pars);
double surface_equator_f(double t, double *xu, void *args);
double surface_equator_df(double t, double *xu, double *xudot, void *args);
double surface_sky_f(double t, double *xu, void *args);
double surface_sky_df(double t, double *xu, double *xudot, void *args);

#endif
