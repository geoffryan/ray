# Ray

A geodesic ray-tracer.

Solves the geodesic equations to trace null worldlines (light rays) from an observer to a target surface through an arbitrary spacetime.

## Dependencies

 - HDF5
 - Lapack

  On Linux, something like the following should get you the dependencies:
```bash
$ sudo apt install libpng-dev lapack-dev libhdf5-dev
```

  On Mac, Lapack is included in the Accelerate framework and should not need to be separately installed, unless you like that sort of thing.  HDF5 can be obtained through Homebrew:
```bash
$ brew install hdf5
```


## Installation

Currently only works out-of-the-box on Mac. Linux coming soon.

Copy `Makefile.in.template` to `Makefile.in` and modify the variables to suit your system (typically just the `H55` variable which should point to your HDF5 installation).  Then just run `make` and you're good to go!

```bash
$ cp Makefile.in.template Makefile.in
$ (... edit Makefile.in...)
$ make
```

## Running

The executable `ray` lives in `bin/`.  It requires a parameter file (such as the example `in.par`) to be provided on the command line at runtime. The parameter file selects the metric tensor, number of rays, position and orientation of the observer, and other run time parameters.

```bash
$ bin/ray in.par
```

This will produce an HDF5 file `map.h5` which contains the skymap of the specified observer, and several `track_????.txt` files which contain the full trajectory data for several individual rays.

## Visualizing

An example visualization script `plotAll.py` lives in `vis/`.  It will produce a plot (`plot.png`) showing projections of the tracks and the observer's view. The default image is a winking smiley face, whose orientation can be used for debugging purposes.
