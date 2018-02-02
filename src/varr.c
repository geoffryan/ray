#include <stdlib.h>
#include "varr.h"

void varr_init(struct varr *v, int size)
{
    v->arr = (double *)malloc(size * sizeof(double));
    v->size = size;
    v->n = 0;
}

void varr_clear(struct varr *v)
{
    v->n = 0;
}

void varr_free(struct varr *v)
{
    if(v->arr != NULL)
    {
        free(v->arr);
        v->arr = NULL;
    }
    v->size = 0;
    v->n = 0;
}

void varr_append(struct varr *v, double x)
{
    int n = v->n;
    if(n >= v->size)
        varr_resize(v, VARR_EXPAND_FAC * v->size);

    v->arr[n] = x;
    (v->n)++;
}

void varr_append_chunk(struct varr *v, double *y, int num)
{
    int n = v->n;
    if(n+num > v->size)
    {
        int new_size = VARR_EXPAND_FAC * v->size;
        while(new_size < n + num)
            new_size *= VARR_EXPAND_FAC;
        varr_resize(v, new_size);
    }

    int i;
    for(i=0; i<num; i++)
        v->arr[n+i] = y[i];
    v->n += num;
}

void varr_resize(struct varr *v, int new_size)
{
    v->arr = (double *)realloc(v->arr, new_size * sizeof(double));
    v->size = new_size;
}

void varr_consolidate(struct varr *v)
{
    varr_resize(v, v->n);
}
