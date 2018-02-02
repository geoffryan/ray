#ifndef RAY_IMAGE
#define RAY_IMAGE

#include "camera.h"

void imageEul(struct Camera *cam, double *x, double tMAX, int nhits, 
                int ntracks, void *args);

#endif
