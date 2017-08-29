#include <stdio.h>
#include <math.h>
#include "par.h"
#include "camera.h"
#include "ode.h"
#include "metric.h"
#include "image.h"
    
void print_tetrad_check(double *e, double *eta, double *g, double *ig);

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        printf("\nusage: geobygeo parfile.par\n");
        printf("exitting\n\n");
        return 0;
    }

    struct parList pars = PAR_DEFAULT;
    struct Camera cam = CAMERA_DEFAULT;
    read_pars(&pars, argv[1]);
    print_pars(&pars, NULL);

    int err = 0;
    err += setup_metric(&pars);
    err += setup_camera(&cam, &pars);

    if(err != 0)
    {
        printf("\nbad setup\n");
        printf("exitting\n\n");
        free_camera(&cam);
        return 0;
    }

    FILE *f = fopen("disc_im_grid.txt", "w");
    fclose(f);
    print_pars(&pars, "disc_im_grid.txt");

    //double args[1];
    //args[0] = 1.0;
    //imageDiscGrid(&grid, args, "disc_im_grid.txt");
    //
    double x[4] = {0.0, 0.0, 0.0, 0.0};
    metric_genPosFromSph(pars.distance, pars.inclination, pars.azimuth, x);
    double args[1] = {1.0};
    double tMAX = 2*pars.distance;

    imageEul(&cam, x, tMAX, args);

    /*
    double g[16], ig[16], eta[16];
    metric_g(g, x, args);
    metric_ig(ig, x, g, args);

    int i,j;
    for(i=0; i<16; i++)
        eta[i] = 0;
    eta[0] = -1;
    eta[5] = 1;
    eta[10] = 1;
    eta[15] = 1;

    double v[3] = {0.1, 0.1, 0.2};
    double u[4];
    u[0] = 1.0/sqrt(-g[0]-2*g[1]*v[0]-2*g[2]*v[1]-2*g[3]*v[2]
                    -g[5]*v[0]*v[0]-2*g[6]*v[0]*v[1]-2*g[7]*v[0]*v[2]
                    -g[10]*v[1]*v[1]-2*g[11]*v[1]*v[2]-g[15]*v[2]*v[2]);
    u[1] = u[0]*v[0];
    u[2] = u[0]*v[1];
    u[3] = u[0]*v[2];

    double norm = 0.0;
    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
            norm += g[4*i+j]*u[i]*u[j];
    printf("u^2: %.6lf\n", norm);
    
    double e[16];

    metric_tetrad_euler(e, x, g, args);

    printf("\nTetrad euler\n\n");
    print_tetrad_check(e, eta, g, ig);

    metric_tetrad(e, x, u, g, args);
    
    printf("\n Tetrad U\n\n");
    print_tetrad_check(e, eta, g, ig);
    */

    free_camera(&cam);

    return 0;
}

void print_tetrad_check(double *e, double *eta, double *g, double *ig)
{
    int i, j, mu, nu;

    printf("e:\n");
    for(i=0; i<4; i++)
    {
        printf("(");
        for(j=0; j<4; j++)
            printf(" %.6lf", e[4*i+j]);
        printf(" )\n");
    }
    printf("ig:\n");
    for(mu=0; mu<4; mu++)
    {
        printf("(");
        for(nu=0; nu<4; nu++)
            printf(" %.6lf", ig[4*mu+nu]);
        printf(" )\n");
    }
    printf("eta^ij e^mu_i e^nu_j:\n");
    for(mu=0; mu<4; mu++)
    {
        printf("(");
        for(nu=0; nu<4; nu++)
        {
            double val = 0.0;
            for(i=0; i<4; i++)
                for(j=0; j<4; j++)
                    val += eta[4*i+j]*e[4*mu+i]*e[4*nu+j];
            printf(" %.6lf", val);
        }
        printf(" )\n");
    }
    printf("eta:\n");
    for(i=0; i<4; i++)
    {
        printf("(");
        for(j=0; j<4; j++)
            printf(" %.6lf", eta[4*i+j]);
        printf(" )\n");
    }
    printf("g_munu e^mu_i e^nu_j:\n");
    for(i=0; i<4; i++)
    {
        printf("(");
        for(j=0; j<4; j++)
        {
            double val = 0.0;
            for(mu=0; mu<4; mu++)
                for(nu=0; nu<4; nu++)
                    val += g[4*mu+nu]*e[4*mu+i]*e[4*nu+j];
            printf(" %.6lf", val);
        }
        printf(" )\n");
    }
}
