#include <stdio.h>
#include "output.h"

void output_print_init(double *args, int n, char *filename)
{
    FILE *f = fopen(filename, "w");
    int i;
    
    for(i=0; i<n; i++)
        fprintf(f, "%.12lg ", args[i]);
    fprintf(f, "\n");

    fclose(f);
}

void output_print_step(double t, double *x, int n, char *filename)
{
    FILE *f = fopen(filename, "a");
    int i;
    fprintf(f, "%.12lg", t);
    for(i=0; i<n; i++)
        fprintf(f, " %.12lg", x[i]);
    fprintf(f, "\n");

    fclose(f);
}

void output_map_h5(double *map, int Nrays, char *filename)
{
    hsize_t fdims2[2];
    hsize_t mapdims[2] = {Nrays, 20};
    hsize_t map_start[2];
    hsize_t f_start[2] = {0,0};

    createFile(filename);
    createGroup(filename, "Rays");

    map_start[0] = 0;

    fdims2[0] = Nrays;
    fdims2[1] = 2;
    createDataset(filename, "Rays", "t", 2, fdims2, H5T_NATIVE_DOUBLE);
    createDataset(filename, "Rays", "thC", 2, fdims2, H5T_NATIVE_DOUBLE);
    map_start[1] = 0;
    writePatch(filename, "Rays", "thC", map, H5T_NATIVE_DOUBLE, 2, 
                map_start, f_start, fdims2, mapdims, fdims2);
    map_start[1] = 2;
    writePatch(filename, "Rays", "t", map, H5T_NATIVE_DOUBLE, 2, 
                map_start, f_start, fdims2, mapdims, fdims2);
    
    fdims2[0] = Nrays;
    fdims2[1] = 4;
    createDataset(filename, "Rays", "x0", 2, fdims2, H5T_NATIVE_DOUBLE);
    createDataset(filename, "Rays", "u0", 2, fdims2, H5T_NATIVE_DOUBLE);
    createDataset(filename, "Rays", "x1", 2, fdims2, H5T_NATIVE_DOUBLE);
    createDataset(filename, "Rays", "u1", 2, fdims2, H5T_NATIVE_DOUBLE);
    map_start[1] = 4;
    writePatch(filename, "Rays", "x0", map, H5T_NATIVE_DOUBLE, 2,
                map_start, f_start, fdims2, mapdims, fdims2);
    map_start[1] = 8;
    writePatch(filename, "Rays", "u0", map, H5T_NATIVE_DOUBLE, 2,
                map_start, f_start, fdims2, mapdims, fdims2);
    map_start[1] = 12;
    writePatch(filename, "Rays", "x1", map, H5T_NATIVE_DOUBLE, 2,
                map_start, f_start, fdims2, mapdims, fdims2);
    map_start[1] = 16;
    writePatch(filename, "Rays", "u1", map, H5T_NATIVE_DOUBLE, 2,
                map_start, f_start, fdims2, mapdims, fdims2);
}

void createFile(char *fname)
{
    hid_t h5file = H5Fcreate( fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5Fclose( h5file );
}

void createGroup(char *fname, char *gname)
{
    hid_t h5file = H5Fopen( fname, H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t h5group = H5Gcreate1( h5file, gname, 0);
    H5Gclose(h5group);
    H5Fclose(h5file);
}

void createDataset(char *fname, char *gname, char *dname, int dim,
                    hsize_t *fdims, hid_t type)
{
    hid_t h5file  = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t h5group = H5Gopen1(h5file, gname);
    hid_t fspace  = H5Screate_simple(dim, fdims, NULL);
    hid_t h5dset  = H5Dcreate1(h5group, dname, type, fspace, H5P_DEFAULT);
    H5Sclose(fspace);
    H5Dclose(h5dset);
    H5Gclose(h5group);
    H5Fclose(h5file);
}

void writeSimple(char *file, char *group, char *dset, void *data, hid_t type)
{
    hid_t h5fil = H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t h5grp = H5Gopen1(h5fil, group);
    hid_t h5dst = H5Dopen1(h5grp, dset);

    H5Dwrite(h5dst, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose(h5dst);
    H5Gclose(h5grp);
    H5Fclose(h5fil);
}

void writePatch(char *file, char *group, char *dset, void *data, hid_t type,
                int dim, hsize_t *loc_start, hsize_t *glo_start, 
                hsize_t *read_size, hsize_t *loc_size, hsize_t *glo_size)
{
    hid_t h5fil = H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t h5grp = H5Gopen1(h5fil, group);
    hid_t h5dst = H5Dopen1(h5grp, dset);

    hsize_t mdims[dim];
    hsize_t fdims[dim];

    hsize_t fstart[dim];
    hsize_t fstride[dim];
    hsize_t fcount[dim];
    hsize_t fblock[dim];
    hsize_t mstart[dim];
    hsize_t mstride[dim];
    hsize_t mcount[dim];
    hsize_t mblock[dim];

    int d;
    for(d=0; d<dim; ++d)
    {
        mdims[d] = loc_size[d];
        fdims[d] = glo_size[d];

        fstart[d]  = glo_start[d];
        fstride[d] = 1;
        fcount[d]  = read_size[d];
        fblock[d]  = 1;
        
        mstart[d]  = loc_start[d];
        mstride[d] = 1;
        mcount[d]  = read_size[d];
        mblock[d]  = 1;
    }
    hid_t mspace = H5Screate_simple(dim, mdims, NULL);
    hid_t fspace = H5Screate_simple(dim, fdims, NULL);

    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, fstart, fstride, fcount, 
                        fblock);
    H5Sselect_hyperslab(mspace, H5S_SELECT_SET, mstart, mstride, mcount, 
                        mblock);

    H5Dwrite(h5dst, type, mspace, fspace, H5P_DEFAULT, data);

    H5Sclose(mspace);
    H5Sclose(fspace);
    H5Dclose(h5dst);
    H5Gclose(h5grp);
    H5Fclose(h5fil);
}
