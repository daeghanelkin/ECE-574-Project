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
