#include <stdio.h>
#include <math.h>
#include "par.h"
#include "camera.h"
#include "output.h"
#include "ode.h"
#include "metric.h"
#include "image.h"
#include "surface.h"
    
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

    output_h5_init(pars.filename);
    output_par_h5(&pars, pars.filename);

    int err = 0;
    err += setup_metric(&pars);
    err += setup_surface(&pars);
    err += setup_camera(&cam, &pars);

    if(err != 0)
    {
        printf("\nbad setup\n");
        printf("exitting\n\n");
        free_camera(&cam);
        return 0;
    }

    double x[4] = {0.0, 0.0, 0.0, 0.0};
    metric_genPosFromSph(pars.distance, pars.inclination, pars.azimuth, x);
    double args[4] = {pars.metricArg1, pars.metricArg2, pars.metricArg3, 
                        pars.metricArg4};
    double tMAX = 2*pars.distance;

    if(pars.observer == 0)
        imageEul(&cam, x, tMAX, pars.nhits, pars.ntracks, args, pars.filename,
                    pars.fancyPrinting);
    else if(pars.observer == 1)
        imageRest(&cam, x, tMAX, pars.nhits, pars.ntracks, args,
                    pars.filename, pars.fancyPrinting);
    else
        printf("\n Bad observer choice: %d", pars.observer);

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
