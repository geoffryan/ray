#include <stdio.h>
#include <math.h>
#ifdef OSX
#include <Accelerate/Accelerate.h>
#endif
#include "par.h"
#include "metric.h"

void metric_print_tensor2(double *t, char name[]);

int setup_metric(struct parList *pars)
{
    int choice = pars->metric;
    int err = 0;

    if(choice == 0)
    {
        metric_orientation = &metric_orientation_flat_cart;
        metric_g  = &metric_g_flat_cart;
        metric_dg = &metric_dg_flat_cart;
        metric_shadow = &metric_shadow_flat_cart;
        metric_fix_domain = &metric_fix_domain_flat_cart;
        metric_genPosFromSph = &metric_genPosFromSph_cart;
    }
    else if(choice == 1)
    {
        metric_orientation = &metric_orientation_flat_sph;
        metric_g  = &metric_g_flat_sph;
        metric_dg = &metric_dg_flat_sph;
        metric_shadow = &metric_shadow_flat_sph;
        metric_fix_domain = &metric_fix_domain_flat_sph;
        metric_genPosFromSph = &metric_genPosFromSph_sph;
    }
    else if(choice == 2)
    {
        metric_orientation = &metric_orientation_schw_sc_sph;
        metric_g  = &metric_g_schw_sc_sph;
        metric_dg = &metric_dg_schw_sc_sph;
        metric_shadow = &metric_shadow_schw_sc_sph;
        metric_fix_domain = &metric_fix_domain_schw_sc_sph;
        metric_genPosFromSph = &metric_genPosFromSph_sph;
    }
    else if(choice == 3)
    {
        metric_orientation = &metric_orientation_schw_ks_sph;
        metric_g  = &metric_g_schw_ks_sph;
        metric_dg = &metric_dg_schw_ks_sph;
        metric_shadow = &metric_shadow_schw_ks_sph;
        metric_fix_domain = &metric_fix_domain_schw_ks_sph;
        metric_genPosFromSph = &metric_genPosFromSph_sph;
    }
    else if(choice == 4)
    {
        metric_orientation = &metric_orientation_schw_ks_cart;
        metric_g  = &metric_g_schw_ks_cart;
        metric_dg = &metric_dg_schw_ks_cart;
        metric_shadow = &metric_shadow_schw_ks_cart;
        metric_fix_domain = &metric_fix_domain_schw_ks_cart;
        metric_genPosFromSph = &metric_genPosFromSph_cart;
    }
    else if(choice == 5)
    {
        metric_orientation = &metric_orientation_bin_ks_cart;
        metric_g  = &metric_g_bin_ks_cart;
        metric_dg = &metric_dg_bin_ks_cart;
        metric_shadow = &metric_shadow_bin_ks_cart;
        metric_fix_domain = &metric_fix_domain_bin_ks_cart;
        metric_genPosFromSph = &metric_genPosFromSph_cart;
    }
    else
    {
        printf("Bad Metric choice: %d\n", choice);
        err = 1;
    }

    return err;
}

void metric_ig(double *ig, double *x, double *g, void *args)
{
    int ipiv[4];
    int n=4;
    int lda=4;
    int info;
    char uplo = 'U';
    int lwork = 256;
    double work[256];

    int i, j;
    for(i=0; i<16; i++)
        ig[i] = g[i];


    //metric_print_tensor2(ig, "g");
    dsytrf_(&uplo, &n, ig, &lda, ipiv, work, &lwork, &info);
    //metric_print_tensor2(ig, "ig_temp");
    dsytri_(&uplo, &n, ig, &lda, ipiv, work, &info);
    //metric_print_tensor2(ig, "ig");

    for(i=1; i<4; i++)
        for(j=2; j<4; j++)
            ig[4*i+j] = ig[4*j+i];
}

void metric_print_tensor2(double *t, char name[])
{
    int i, j;
    printf("%s\n", name);
    for(i=0; i<4; i++)
    {
        printf("( ");
        for(j=0; j<4; j++)
            printf("%.8lf ", t[4*i+j]);
        printf(")\n");
    }
}

void metric_genPosFromSph_cart(double dist, double inc, double az, double *X)
{
    double x = dist*sin(M_PI*inc/180.0)*cos(M_PI*az/180.0);
    double y = dist*sin(M_PI*inc/180.0)*sin(M_PI*az/180.0);
    double z = dist*cos(M_PI*inc/180.0);

    X[1] = x;
    X[2] = y;
    X[3] = z;
}

void metric_genPosFromSph_sph(double dist, double inc, double az, double *X)
{
    X[1] = dist;
    X[2] = M_PI*inc/180.0;
    X[3] = M_PI*az/180.0;
}
void metric_tetrad_euler(double *e, double *x, double *g, 
                            void *args)
{
    double ig[16];
    metric_ig(ig, x, g, args);


    double e0[4], e1[4], e2[4], e3[4];

    double al = sqrt(-1.0/ig[0]);
    double be[3] = {-ig[1]/ig[0], -ig[2]/ig[0], -ig[3]/ig[0]};
    e0[0] = 1.0/al;
    e0[1] = -be[0]/al;
    e0[2] = -be[1]/al;
    e0[3] = -be[2]/al;

    int mu;

    for(mu=0; mu<4; mu++)
    {
        e1[mu] = 0;
        e2[mu] = 0;
        e3[mu] = 0;
    }

    // e3 parallel to x3
    e3[3] = 1.0/sqrt(g[15]);

    // Gram-Schmidt for e2
    e2[2] = 1.0;
    e2[3] = 1.0;
    double b2de3 = e3[3]*(g[14]*e2[2]+g[15]*e2[3]);
    e2[3] -= b2de3*e3[3];
    double norm = sqrt(g[10]*e2[2]*e2[2] + 2*g[11]*e2[2]*e2[3]
                        + g[15]*e2[3]*e2[3]);
    e2[2] /= norm;
    e2[3] /= norm;
    
    // Gram-Schmidt for e1
    e1[1] = 1.0;
    e1[2] = 1.0;
    e1[3] = 1.0;
    double b1de3 = e3[3]*(g[13]*e1[1] + g[14]*e1[2] + g[15]*e1[3]);
    double b1de2 = e2[2]*(g[9]*e1[1] + g[10]*e1[2] + g[11]*e1[3])
                    + e2[3]*(g[13]*e1[1] + g[14]*e1[2] + g[15]*e1[3]);
    e1[2] -= b1de2*e2[2];
    e1[3] -= b1de2*e2[3] + b1de3*e3[3];
    norm = sqrt(g[5]*e1[1]*e1[1] + 2*g[6]*e1[1]*e1[2]
                        + 2*g[6]*e1[1]*e1[3] + g[10]*e1[2]*e1[2]
                        + 2*g[11]*e1[2]*e1[3] + g[15]*e1[3]*e1[3]);
    e1[1] /= norm;
    e1[2] /= norm;
    e1[3] /= norm;

    for(mu=0; mu<4; mu++)
    {
        e[4*mu+0] = e0[mu];
        e[4*mu+1] = e1[mu];
        e[4*mu+2] = e2[mu];
        e[4*mu+3] = e3[mu];
    }
}

void metric_tetrad(double *e, double *x, double *u, double *g, void *args)
{
    double ee[16];
    metric_tetrad_euler(ee, x, g, args);

    double iee[16];
    int i, j, mu, nu;
    double eta[16];
    for(i=0; i<16; i++)
        eta[i] = 0.0;
    eta[0] = -1;
    eta[5] = 1;
    eta[10] = 1;
    eta[15] = 1;

    for(i=0; i<4; i++)
        for(mu=0; mu<4; mu++)
        {
            iee[4*i+mu] = 0.0;
            for(nu=0; nu<4; nu++)
                for(j=0; j<4; j++)
                    iee[4*i+mu] += g[4*mu+nu]*eta[4*i+j]*ee[4*nu+j];
        }

    double ue[4] = {0., 0., 0., 0.};
    for(i=0; i<4; i++)
        for(mu=0; mu<4; mu++)
            ue[i] += iee[4*i+mu]*u[mu];

    //printf("ue: %.6lf %.6lf %.6lf %.6lf\n", ue[0], ue[1], ue[2], ue[3]);

    double gam = ue[0];
    double u2 = ue[1]*ue[1] + ue[2]*ue[2] + ue[3]*ue[3];
    double v2 = u2 / (gam*gam);
    double v = sqrt(v2);
    double gm1 = u2/(gam+1);

    double n1 = 1.0;
    double n2 = 0.0;
    double n3 = 0.0;

    if(v > 1.0e-15)
    {
        n1 = ue[1] / (gam*v);
        n2 = ue[2] / (gam*v);
        n3 = ue[3] / (gam*v);
    }

    //printf("gam: %.6lf u2: %.6lf v2: %.6lf v: %.6lf gam-1: %.6lf\n",
    //        gam, u2, v2, v, gm1);
    //printf("n: %.6lf %.6lf %.6lf\n", n1, n2, n3);
    
    //Lorentz matrix to boost into u's frame.
    double lor[16];
    lor[4*0+0] = gam;
    lor[4*0+1] = -gam*v*n1;
    lor[4*0+2] = -gam*v*n2;
    lor[4*0+3] = -gam*v*n3;
    lor[4*1+0] = lor[4*0+1];
    lor[4*1+1] = 1 + gm1*n1*n1;
    lor[4*1+2] = gm1*n1*n2;
    lor[4*1+3] = gm1*n1*n3;
    lor[4*2+0] = lor[4*0+2];
    lor[4*2+1] = lor[4*1+2];
    lor[4*2+2] = 1 + gm1*n2*n2;
    lor[4*2+3] = gm1*n2*n3;
    lor[4*3+0] = lor[4*0+3];
    lor[4*3+1] = lor[4*1+3];
    lor[4*3+2] = lor[4*2+3];
    lor[4*3+3] = 1 + gm1*n3*n3;

    //dual Lorentz matrix to transform ee
    lor[1] = -lor[1];
    lor[2] = -lor[2];
    lor[3] = -lor[3];
    lor[4] = -lor[4];
    lor[8] = -lor[8];
    lor[12] = -lor[12];

    for(mu=0; mu<4; mu++)
        for(i=0; i<4; i++)
        {
            e[4*mu+i] = 0.0;
            for(j=0; j<4; j++)
                e[4*mu+i] += lor[4*i+j]*ee[4*mu+j];
        }

    //printf("u:  %.10lg %.10lg %.10lg %.10lg\n", u[0], u[1], u[2], u[3]);
    //printf("e0: %.10lg %.10lg %.10lg %.10lg\n", e[0], e[4], e[8], e[12]);
}

