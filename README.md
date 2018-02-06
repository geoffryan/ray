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

This will produce an HDF5 file (by default `map.h5`) which contains the skymap of the specified observer, a copy of the runtime paramerters, and several tracks which contain the full trajectory data for several individual rays.

## Visualizing

Utilities for interacting with and visualizing the `ray` data live in `vis/`, including example visualization scripts `plotAll.py`, `plotFace.py`, and `plotMovingFace.py`.  All can be run on the command line and take as input a `ray` HDF5 file. 

```bash
$ python vis/plotAll.py map.h5 track_*.txt
```

`plotAll.py` produces a plot (`plot.png`) showing projections of the tracks and the observer's view. The default image is a winking smiley face, whose orientation can be used for debugging purposes (open eye is -x,+y quadrant, closed eye +x,+y quadrant, and the smile on the -y half plane).

`plotFace.py` and `plotMovingFace.py` only produce observer images of the smiley face.


## Acknowledgements

Thanks to Micha Gorelick for testing on Linux and making me add this README.
