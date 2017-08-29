#include "metric.h"
#include "trace.h"
#include "image.h"

void imageEul(struct Camera *cam, double *x, double tMAX, void *args)
{
    double g[16], e[16];
    metric_g(g, x, args);
    metric_tetrad_euler(e, x, g, args);
    double u[] = {e[0], e[4], e[8], e[12]};

    camera_set_position(cam, x, u, args);

    trace(cam, tMAX, 1617, target_eq_sph, args);
}
