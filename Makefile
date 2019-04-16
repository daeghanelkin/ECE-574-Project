CC = gcc
MCC = mpicc
CFLAGS = -O2 -g -Wall
LFLAGS = -lX11 -ljpeg -lm

all: sobel_before sobel_mpi

sobel_before: sobel_before.o
	$(CC) -o sobel_before sobel_before.o $(LFLAGS)

sobel_before.o: sobel_before.c
	$(CC) $(CFLAGS) -c sobel_before.c 

sobel_mpi: sobel_mpi.o
	$(MCC) -o sobel_mpi sobel_mpi.o $(LFLAGS)

sobel_mpi.o: sobel_mpi.c
	$(MCC) $(CFLAGS) -c sobel_mpi.c

clean:
	rm -f *~ *.o sobel_before sobel_mpi out.jpg
