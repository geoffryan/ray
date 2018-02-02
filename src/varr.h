#ifndef RAY_VARARR
#define RAY_VARARR

#define VARR_EXPAND_FAC 1.5

struct varr
{
    int n;
    int size;
    double *arr;
};

const static struct varr VARR_DEFAULT = {0, 0, NULL};

void varr_init(struct varr *v, int size);
void varr_clear(struct varr *v);
void varr_free(struct varr *v);
void varr_append(struct varr *v, double x);
void varr_append_chunk(struct varr *v, double *y, int num);
void varr_resize(struct varr *v, int new_size);
void varr_consolidate(struct varr *v);

#endif
