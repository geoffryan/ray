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

void output_h5_init(char *filename)
{
    createFile(filename);
    createGroup(filename, "Map");
    createGroup(filename, "Tracks");
    createGroup(filename, "Pars");
}

void output_par_h5(struct parList *pars, char *filename)
{
    hsize_t fdims1[1] = {1};
    
    createDataset(filename, "Pars", "Metric", 1, fdims1, H5T_NATIVE_INT);
    writeSimple(filename, "Pars", "Metric", &(pars->metric), H5T_NATIVE_INT);
    
    createDataset(filename, "Pars", "MeshType", 1, fdims1, H5T_NATIVE_INT);
    writeSimple(filename, "Pars", "MeshType", &(pars->meshType),
                                                            H5T_NATIVE_INT);
    createDataset(filename, "Pars", "N1", 1, fdims1, H5T_NATIVE_INT);
    writeSimple(filename, "Pars", "N1", &(pars->N1), H5T_NATIVE_INT); 
    createDataset(filename, "Pars", "N2", 1, fdims1, H5T_NATIVE_INT);
    writeSimple(filename, "Pars", "N2", &(pars->N2), H5T_NATIVE_INT);
    
    createDataset(filename, "Pars", "X1a", 1, fdims1, H5T_NATIVE_DOUBLE);
    writeSimple(filename, "Pars", "X1a", &(pars->X1a), H5T_NATIVE_DOUBLE); 
    createDataset(filename, "Pars", "X1b", 1, fdims1, H5T_NATIVE_DOUBLE);
    writeSimple(filename, "Pars", "X1b", &(pars->X1b), H5T_NATIVE_DOUBLE); 
    createDataset(filename, "Pars", "X2a", 1, fdims1, H5T_NATIVE_DOUBLE);
    writeSimple(filename, "Pars", "X2a", &(pars->X2a), H5T_NATIVE_DOUBLE); 
    createDataset(filename, "Pars", "X2b", 1, fdims1, H5T_NATIVE_DOUBLE);
    writeSimple(filename, "Pars", "X2b", &(pars->X2b), H5T_NATIVE_DOUBLE); 

    createDataset(filename, "Pars", "nhits", 1, fdims1, H5T_NATIVE_INT);
    writeSimple(filename, "Pars", "nhits", &(pars->nhits), H5T_NATIVE_INT); 

    createDataset(filename, "Pars", "distance", 1, fdims1, H5T_NATIVE_DOUBLE);
    writeSimple(filename, "Pars", "distance", &(pars->distance), 
                                                        H5T_NATIVE_DOUBLE); 
    createDataset(filename, "Pars", "inclination", 1, fdims1,
                                                        H5T_NATIVE_DOUBLE);
    writeSimple(filename, "Pars", "inclination", &(pars->inclination), 
                                                        H5T_NATIVE_DOUBLE); 
    createDataset(filename, "Pars", "azimuth", 1, fdims1, H5T_NATIVE_DOUBLE);
    writeSimple(filename, "Pars", "azimuth", &(pars->azimuth), 
                                                        H5T_NATIVE_DOUBLE); 
}

void output_map_h5(double *map, int Nrays, char *filename)
{
    hsize_t fdims2[2];
    hsize_t mapdims[2] = {Nrays, 20};
    hsize_t map_start[2];
    hsize_t f_start[2] = {0,0};

    map_start[0] = 0;

    fdims2[0] = Nrays;
    fdims2[1] = 2;
    createDataset(filename, "Map", "t", 2, fdims2, H5T_NATIVE_DOUBLE);
    createDataset(filename, "Map", "thC", 2, fdims2, H5T_NATIVE_DOUBLE);
    map_start[1] = 0;
    writePatch(filename, "Map", "thC", map, H5T_NATIVE_DOUBLE, 2, 
                map_start, f_start, fdims2, mapdims, fdims2);
    map_start[1] = 2;
    writePatch(filename, "Map", "t", map, H5T_NATIVE_DOUBLE, 2, 
                map_start, f_start, fdims2, mapdims, fdims2);
    
    fdims2[0] = Nrays;
    fdims2[1] = 4;
    createDataset(filename, "Map", "x0", 2, fdims2, H5T_NATIVE_DOUBLE);
    createDataset(filename, "Map", "u0", 2, fdims2, H5T_NATIVE_DOUBLE);
    createDataset(filename, "Map", "x1", 2, fdims2, H5T_NATIVE_DOUBLE);
    createDataset(filename, "Map", "u1", 2, fdims2, H5T_NATIVE_DOUBLE);
    map_start[1] = 4;
    writePatch(filename, "Map", "x0", map, H5T_NATIVE_DOUBLE, 2,
                map_start, f_start, fdims2, mapdims, fdims2);
    map_start[1] = 8;
    writePatch(filename, "Map", "u0", map, H5T_NATIVE_DOUBLE, 2,
                map_start, f_start, fdims2, mapdims, fdims2);
    map_start[1] = 12;
    writePatch(filename, "Map", "x1", map, H5T_NATIVE_DOUBLE, 2,
                map_start, f_start, fdims2, mapdims, fdims2);
    map_start[1] = 16;
    writePatch(filename, "Map", "u1", map, H5T_NATIVE_DOUBLE, 2,
                map_start, f_start, fdims2, mapdims, fdims2);
}

void output_track_h5(struct varr *track, int id, char *filename)
{
    char gname[128];
    sprintf(gname, "Tracks/track_%06d", id);

    int nstep = track->n / 9;
    hsize_t fdims1[1] = {nstep};
    hsize_t fdims2[2] = {nstep,4};

    createGroup(filename, gname);
    createDataset(filename, gname, "t", 1, fdims1, H5T_NATIVE_DOUBLE);
    createDataset(filename, gname, "x", 2, fdims2, H5T_NATIVE_DOUBLE);
    createDataset(filename, gname, "u", 2, fdims2, H5T_NATIVE_DOUBLE);
   
    hsize_t fstart2[2] = {0,0};
    hsize_t track_start2[2] = {0,0};
    hsize_t track_dims2[2] = {nstep,9};

    fdims2[1] = 1;
    track_start2[1] = 0;
    writePatch(filename, gname, "t", track->arr, H5T_NATIVE_DOUBLE, 2,
                track_start2, fstart2, fdims2, track_dims2, fdims2);

    fdims2[1] = 4;
    track_start2[1] = 1;
    writePatch(filename, gname, "x", track->arr, H5T_NATIVE_DOUBLE, 2,
                track_start2, fstart2, fdims2, track_dims2, fdims2);
    track_start2[1] = 5;
    writePatch(filename, gname, "u", track->arr, H5T_NATIVE_DOUBLE, 2,
                track_start2, fstart2, fdims2, track_dims2, fdims2);
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
