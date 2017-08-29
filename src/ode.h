#ifndef RAY_ODE
#define RAY_ODE

double forward_euler(double t, double *x, double *dx, int n, void *args,
                    void (*xdot)(double,double*,void*,double*));
double rk2(double t, double *x, double *dx, int n, void *args,
                    void (*xdot)(double,double*,void*,double*));
double rk4(double t, double *x, double *dx, int n, void *args,
                    void (*xdot)(double,double*,void*,double*));
double dopr54(double t, double *x, double *dt, int n, void *args,
                    void (*xdot)(double,double*,void*,double*));

#endif
