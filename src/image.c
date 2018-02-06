#include <math.h>
#include <stdio.h>
#include "metric.h"
#include "trace.h"
#include "image.h"

void imageEul(struct Camera *cam, double *x, double tMAX, int nhits,
                int ntracks, void *args, char *filename)
{
    double g[16], e[16];
    metric_g(g, x, args);
    metric_tetrad_euler(e, x, g, args);
    double u[] = {e[0], e[4], e[8], e[12]};

    camera_set_position(cam, x, u, args);
    camera_print(cam);

    trace(cam, tMAX, nhits, ntracks, args, filename);
}

void imageRest(struct Camera *cam, double *x, double tMAX, int nhits,
                int ntracks, void *args, char *filename)
{
    double g[16];
    metric_g(g, x, args);

    if(g[0] >= 0.0)
    {
        printf("Observer has unphysical velocity\n.");
        return;
    }

    double u[] = {sqrt(-1/g[0]), 0.0, 0.0, 0.0};
    camera_set_position(cam, x, u, args);
    camera_print(cam);

    trace(cam, tMAX, nhits, ntracks, args, filename);
}
