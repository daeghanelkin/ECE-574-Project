CC = gcc
CFLAGS = -g -Wall
LIBS = -lX11 -ljpeg -lm

all: 
	$(CC) -o sobel_before $(CFLAGS) sobel_before.c $(LIBS)

clean:
	rm -f sobel_before
