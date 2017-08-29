#include <math.h>
#include <stdlib.h>
#include "ode.h"

double forward_euler(double t, double *x, double *dt, int n, void *args,
                    void (*xdot)(double,double*,void*,double*))
{
    // Simple forward euler integration scheme.  Updates 'x' in place.

    int i;
    double *k = (double *) malloc(n * sizeof(double));
    
    xdot(t, x, args, k);

    for(i=0; i<n; i++)
        x[i] += (*dt)*k[i];

    free(k);

    return *dt;
}

double rk2(double t, double *x, double *dt, int n, void *args,
            void (*xdot)(double,double*,void*,double*))
{
    // Basic Runge-Kutta 2nd order integration scheme.  Updates 'x' in place.

    int i;
    double *k1 = (double *) malloc(n * sizeof(double));
    double *k2 = (double *) malloc(n * sizeof(double));
    double *x1 = (double *) malloc(n * sizeof(double));
    
    xdot(t, x, args, k1);

    for(i=0; i<n; i++)
        x1[i] = x[i] + 0.5*(*dt)*k1[i];
    xdot(t, x1, args, k2);

    for(i=0; i<n; i++)
        x[i] += (*dt)*k2[i];

    free(k1);
    free(k2);
    free(x1);
    
    return *dt;
}

double rk4(double t, double *x, double *dt, int n, void *args,
            void (*xdot)(double,double*,void*,double*))
{
    // Basic Runge-Kutta 4th order integration scheme.  Updates 'x' in place.

    int i;
    double *k1 = (double *) malloc(n * sizeof(double));
    double *k2 = (double *) malloc(n * sizeof(double));
    double *k3 = (double *) malloc(n * sizeof(double));
    double *k4 = (double *) malloc(n * sizeof(double));
    double *x1 = (double *) malloc(n * sizeof(double));
    double *x2 = (double *) malloc(n * sizeof(double));
    double *x3 = (double *) malloc(n * sizeof(double));
    
    xdot(t, x, args, k1);

    for(i=0; i<n; i++)
        x1[i] = x[i] + 0.5*(*dt)*k1[i];
    xdot(t, x1, args, k2);

    for(i=0; i<n; i++)
        x2[i] = x[i] + 0.5*(*dt)*k2[i];
    xdot(t, x2, args, k3);

    for(i=0; i<n; i++)
        x3[i] = x[i] + (*dt)*k3[i];
    xdot(t, x3, args, k4);

    for(i=0; i<n; i++)
        x[i] += (*dt)*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(x1);
    free(x2);
    free(x3);

    return *dt;
}

double dopr54(double t, double *x, double *dt, int n, void *args,
            void (*xdot)(double,double*,void*,double*))
{
    // Dormand-Prince 4/5'th order RK method with adaptive step sizing.
    // Updates 'x' in place.
        
    double atol = 1.0e-10;
    double rtol = 1.0e-10;
    
    int i;
    double *k1 = (double *) malloc(n * sizeof(double));
    double *k2 = (double *) malloc(n * sizeof(double));
    double *k3 = (double *) malloc(n * sizeof(double));
    double *k4 = (double *) malloc(n * sizeof(double));
    double *k5 = (double *) malloc(n * sizeof(double));
    double *k6 = (double *) malloc(n * sizeof(double));
    double *k7 = (double *) malloc(n * sizeof(double));
    double *x1 = (double *) malloc(n * sizeof(double));
    double *x2 = (double *) malloc(n * sizeof(double));
    double *x3 = (double *) malloc(n * sizeof(double));
    double *x4 = (double *) malloc(n * sizeof(double));
    double *x5 = (double *) malloc(n * sizeof(double));
    double *x6 = (double *) malloc(n * sizeof(double));
    double *xf5 = (double *) malloc(n * sizeof(double));
    double *xf4 = (double *) malloc(n * sizeof(double));
   
    double b1 = 35.0/384.0;
    double b2 = 0.0;
    double b3 = 500.0/1113.0;
    double b4 = 125.0/192.0;
    double b5 = -2187.0/6784.0;
    double b6 = 11.0/84.0;
    double b7 = 0.0;

    double d1 = 5179.0/57600.0;
    double d2 = 0.0;
    double d3 = 7571.0/16695.0;
    double d4 = 393.0/640.0;
    double d5 = -92097.0/339200.0;
    double d6 = 187.0/2100.0;
    double d7 = 1.0/40.0;

    double c2 = 0.2;
    double c3 = 0.3;
    double c4 = 0.8;
    double c5 = 8.0/9.0;
    double c6 = 1.0;
    double c7 = 1.0;

    double a21 = 0.2;
    double a31 = 3.0/40.0;
    double a32 = 9.0/40.0;
    double a41 = 44.0/45.0;
    double a42 = -56.0/15.0;
    double a43 = 32.0/9.0;
    double a51 = 19372.0/6561.0;
    double a52 = -25360.0/2187.0;
    double a53 = 64448.0/6561.0;
    double a54 = -212.0/729.0;
    double a61 = 9017.0/3168.0;
    double a62 = -355.0/33.0;
    double a63 = 46732.0/5247.0;
    double a64 = 49.0/176.0;
    double a65 = -5103.0/18656.0;
    double a71 = 35.0/384.0;
    double a72 = 0.0;
    double a73 = 500.0/1113.0;
    double a74 = 125.0/192.0;
    double a75 = -2187.0/6784.0;
    double a76 = 11.0/84.0;

    double err = 1.0e100;
    double dt0 = *dt;
    double dt1 = *dt;

    xdot(t, x, args, k1);

    while(err > 1.0)
    {
        dt0 = dt1;

        for(i=0; i<n; i++)
            x1[i] = x[i] + dt0*(a21*k1[i]);
        xdot(t+c2*dt0, x1, args, k2);

        for(i=0; i<n; i++)
            x2[i] = x[i] + dt0*(a31*k1[i]+a32*k2[i]);
        xdot(t+c3*dt0, x2, args, k3);

        for(i=0; i<n; i++)
            x3[i] = x[i] + dt0*(a41*k1[i]+a42*k2[i]+a43*k3[i]);
        xdot(t+c4*dt0, x3, args, k4);

        for(i=0; i<n; i++)
            x4[i] = x[i] + dt0*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]);
        xdot(t+c5*dt0, x4, args, k5);

        for(i=0; i<n; i++)
            x5[i] = x[i] + dt0*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]
                                    +a65*k5[i]);
        xdot(t+c6*dt0, x5, args, k6);

        for(i=0; i<n; i++)
            x6[i] = x[i] + dt0*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]
                                    +a75*k5[i]+a76*k6[i]);
        xdot(t+c7*dt0, x6, args, k7);

        int i;
        for(i=0; i<n; i++)
        {
            xf5[i] = x[i] + dt0*(b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i]
                                    + b5*k5[i] + b6*k6[i] + b7*k7[i]);
            xf4[i] = x[i] + dt0*(d1*k1[i] + d2*k2[i] + d3*k3[i] + d4*k4[i]
                                    + d5*k5[i] + d6*k6[i] + d7*k7[i]);
        }

        err = 0.0;
        double xm;

        for(i=0; i<n; i++)
        {
            xm = fabs(x[i]) > fabs(xf5[i]) ? fabs(x[i]) : fabs(xf5[i]);
            err += (xf5[i]-xf4[i])*(xf5[i]-xf4[i])
                    / ((atol+xm*rtol)*(atol+xm*rtol));
        }
        err = sqrt(err/n);
        if(err > 0.0)
            dt1 = 0.9 * dt0 * pow(err, -0.2);
        else
            dt1 = 2.0*dt0;
    }

    *dt = dt1;

    for(i=0; i<n; i++)
        x[i] = xf5[i];


    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(k5);
    free(k6);
    free(k7);
    free(x1);
    free(x2);
    free(x3);
    free(x4);
    free(x5);
    free(x6);
    free(xf4);
    free(xf5);

    return dt0;
}

