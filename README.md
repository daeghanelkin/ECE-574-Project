# ECE-574-Project

## Overview

This is a project for the University of Maine Electrical and Computer Engineering Cluster Computing Graduate course. 
The goal of this project is to run sobel image processing code on a Raspberry Pi cluster using MPI while providing a
X11 interface to watch the image being processed.

This was intended to be run on Professor Vince Weaver's Pi Cluster to show to anyone touring or just wants to see 
a visualization of what a cluster can do.

## Setup

There are a few required libraries for compiling the program. 
* libjpeg
* libx11 or libx11-dev (Depending on system you're running on)
* mpi or openmpi

MPI might require some setup of its own so head to their [site](https://www.open-mpi.org/) for that.

## Compiling

Simply running `make` will create the executables for both `sobel_before` and `sobel_mpi`.

If you wish to compile a single source simply run `make sobel_before` or `make sobel_mpi`.

Running `make clean` will clean the directory of anything not required for compiling/running.

## Running

Running the executables can be done two ways.

For `sobel_before` you can simply run `./sobel_before <input_image.jpg>`.

For `sobel_mpi` you need to run `mpiexec -n <# of Cores> ./sobel_mpi <input_image.jpg>`.

## Future Work

This project isn't anywhere near perfect. Some future improvements might be with
performance, image sizing, and information. This run very slow, even with small
images so improvements to how we gather or divide up the image per core might
improve that. Currently this has no support for images larger than your display.
If an image is larger than your resolution you'll only see part of it. There's
probably some way to scale but that's for another time. It would also be cool to
show some information about work load/distribution across the pi cluster.

