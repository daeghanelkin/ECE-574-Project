CC = gcc
CFLAGS = -g -Wall
LIBS = -lX11 -ljpeg

all: 
	$(CC) -o sobel_before $(CFLAGS) sobel_before.c $(LIBS)

clean:
	rm -f sobel_before