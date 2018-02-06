#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "par.h"

int readvar(char filename[], char key[], int vtype, void *ptr)
{
    FILE *f = fopen(filename, "r");

    char line[MAXSTRLEN];
    char word[MAXSTRLEN];
    int found = 0;

    while(fgets(line,MAXSTRLEN,f) != NULL)
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

    int whitespace_length = strspn(line+strlen(key)," \t:=");
    char *sval = line + strlen(key) + whitespace_length;

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
        int entry_length = strcspn(sval, "\t:=#\n");
        char entry[MAXSTRLEN];
        strncpy(entry, sval, entry_length);
        entry[entry_length] = '\0';

        char *comment = strstr(entry, "//");
        if(comment != NULL)
            *comment = '\0';
        entry_length = strlen(entry);

        comment = strstr(entry, "/*");
        if(comment != NULL)
            *comment = '\0';
        entry_length = strlen(entry);

        //Trim trailing whitespace
        char *end = entry + entry_length - 1;
        while(end > entry && isspace(*end))
            end--;
        *(end+1) = '\0';
        entry_length = strlen(entry);

        strncpy((char *) ptr, entry, entry_length);
        ((char *)ptr)[entry_length] = '\0';
    }

    return 0;
}

void read_pars(struct parList *theParList, char filename[])
{
    readvar(filename, "Filename",       VAR_STR, &(theParList->filename));
    readvar(filename, "Metric",          VAR_INT, &(theParList->metric));
    readvar(filename, "Surface",          VAR_INT, &(theParList->surface));
    readvar(filename, "Observer",       VAR_INT, &(theParList->observer));
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
    readvar(filename, "FancyPrinting",  VAR_INT, &(theParList->fancyPrinting));
}

void print_pars(struct parList *theParList, char filename[])
{
    FILE *f;
    if(filename == NULL)
        f = stdout;
    else
        f = fopen(filename, "a");

    fprintf(f, "### Input Parameters ###\n");
    fprintf(f, "Filename: %s\n", theParList->filename);
    fprintf(f, "Metric: %d\n", theParList->metric);
    fprintf(f, "Surface: %d\n", theParList->surface);
    fprintf(f, "Observer: %d\n", theParList->observer);
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
    fprintf(f, "fancyPrinting: %d\n", theParList->fancyPrinting);

    if(filename != NULL)
        fclose(f);
}
