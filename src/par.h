#ifndef RAY_PAR
#define RAY_PAR

#define MAXSTRLEN 256

enum{VAR_DBL, VAR_INT, VAR_LON, VAR_STR};

struct parList
{
    char filename[MAXSTRLEN];
    int metric;
    int surface;
    int observer;
    int meshType;
    int N1;
    int N2;
    int nhits;
    int ntracks;
    double metricArg1;
    double metricArg2;
    double metricArg3;
    double metricArg4;
    double X1a;
    double X1b;
    double X2a;
    double X2b;
    double distance;
    double inclination;
    double azimuth;
    int fancyPrinting;
};

const static struct parList PAR_DEFAULT = {
    .filename = "map.h5",
    .metric = 0,
    .surface = 0,
    .observer = 0,
    .meshType = 0,
    .N1 = 1,
    .N2 = 1,
    .nhits = 1,
    .ntracks = 100,
    .metricArg1 = 0.0,
    .metricArg2 = 0.0,
    .metricArg3 = 0.0,
    .metricArg4 = 0.0,
    .X1a = 0.0,
    .X1b = 0.0,
    .X2a = 0.0,
    .X2b = 0.0,
    .distance = 0.0,
    .inclination = 0.0,
    .azimuth = 0.0,
    .fancyPrinting = 0
};

int readvar(char filename[], char key[], int vtype, void *ptr);
void read_pars(struct parList *theParList, char filename[]);
void print_pars(struct parList *theParList, char filename[]);

#endif
