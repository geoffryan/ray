#include <stdio.h>
#include <string.h>
#include "par.h"

int readvar(char filename[], char key[], int vtype, void *ptr)
{
    FILE *f = fopen(filename, "r");

    char line[256];
    char word[256];
    int found = 0;

    while(fgets(line,256,f) != NULL)
    {
        sscanf(line, "%s ", word);
        if(strcmp(word,key) == 0)
        {
            found = 1;
            break;
        }
    }
    fclose(f);
    if(!found)
    {
        printf("SETUP: %s parameter not found. Using default.\n", key);
        return 1;
    }

    char *sval = line + strlen(key) + strspn(line+strlen(key)," \t:=");

    if(vtype == VAR_DBL)
    {
        double val;
        sscanf(sval, "%lf", &val);
        *((double *)ptr) = val;
    }
    else if(vtype == VAR_INT)
    {
        int val;
        sscanf(sval, "%d", &val);
        *((int *)ptr) = val;
    }
    else if(vtype == VAR_LON)
    {
        long val;
        sscanf(sval, "%ld", &val);
        *((long *)ptr) = val;
    }
    else
    {
        strcpy((char *) ptr, sval);
    }

    return 0;
}

void read_pars(struct parList *theParList, char filename[])
{
    readvar(filename, "Metric",          VAR_INT, &(theParList->metric));
    readvar(filename, "MeshType",          VAR_INT, &(theParList->meshType));
    readvar(filename, "N1",          VAR_INT, &(theParList->N1));
    readvar(filename, "N2",          VAR_INT, &(theParList->N2));
    readvar(filename, "nhits",        VAR_INT, &(theParList->nhits));
    readvar(filename, "ntracks",        VAR_INT, &(theParList->ntracks));
    readvar(filename, "X1a",          VAR_DBL, &(theParList->X1a));
    readvar(filename, "X1b",          VAR_DBL, &(theParList->X1b));
    readvar(filename, "X2a",          VAR_DBL, &(theParList->X2a));
    readvar(filename, "X2b",          VAR_DBL, &(theParList->X2b));
    readvar(filename, "inclination",  VAR_DBL, &(theParList->inclination));
    readvar(filename, "azimuth",  VAR_DBL, &(theParList->azimuth));
    readvar(filename, "distance",  VAR_DBL, &(theParList->distance));
    readvar(filename, "imDim",  VAR_INT, &(theParList->imDim));
    readvar(filename, "imX",  VAR_DBL, &(theParList->imX));
    readvar(filename, "imPiFac",  VAR_INT, &(theParList->imPiFac));
}

void print_pars(struct parList *theParList, char filename[])
{
    FILE *f;
    if(filename == NULL)
        f = stdout;
    else
        f = fopen(filename, "a");

    fprintf(f, "### Input Parameters ###\n");
    fprintf(f, "Metric: %d\n", theParList->metric);
    fprintf(f, "MeshType: %d\n", theParList->meshType);
    fprintf(f, "N1: %d\n", theParList->N1);
    fprintf(f, "N2: %d\n", theParList->N2);
    fprintf(f, "nhits: %d\n", theParList->nhits);
    fprintf(f, "ntracks: %d\n", theParList->ntracks);
    fprintf(f, "X1a: %g\n", theParList->X1a);
    fprintf(f, "X1b: %g\n", theParList->X1b);
    fprintf(f, "X2a: %g\n", theParList->X2a);
    fprintf(f, "X2b: %g\n", theParList->X2b);
    fprintf(f, "distance: %g\n", theParList->distance);
    fprintf(f, "inclination: %g\n", theParList->inclination);
    fprintf(f, "azimuth: %g\n", theParList->azimuth);
    fprintf(f, "imDim: %d\n", theParList->imDim);
    fprintf(f, "imX: %g\n", theParList->imX);
    fprintf(f, "imPiFac: %d\n", theParList->imPiFac);

    if(filename != NULL)
        fclose(f);
}
