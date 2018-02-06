#ifndef RAY_IMAGE
#define RAY_IMAGE

#include "camera.h"

void imageEul(struct Camera *cam, double *x, double tMAX, int nhits, 
                int ntracks, void *args, char *filename, int fancyPrinting);
void imageRest(struct Camera *cam, double *x, double tMAX, int nhits,
                int ntracks, void *args, char *filename, int fancyPrinting);

#endif
