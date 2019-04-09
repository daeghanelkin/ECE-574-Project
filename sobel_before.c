/* Example sobel code for ECE574 -- Spring 2017 */
/* By Vince Weaver <vincent.weaver@maine.edu> */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>

#include <jpeglib.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

/* Structure describing the image */
struct image_t {
	int x;
	int y;
	int depth;	/* bytes */
	unsigned char *pixels;
};

static int load_jpeg(char *filename, struct image_t *image) {

	FILE *fff;
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPROW output_data;
	unsigned int scanline_len;
	int scanline_count=0;

	fff=fopen(filename,"rb");
	if (fff==NULL) {
		fprintf(stderr, "Could not load %s: %s\n",
			filename, strerror(errno));
		return -1;
	}

	/* set up jpeg error routines */
	cinfo.err = jpeg_std_error(&jerr);

	/* Initialize cinfo */
	jpeg_create_decompress(&cinfo);

	/* Set input file */
	jpeg_stdio_src(&cinfo, fff);

	/* read header */
	jpeg_read_header(&cinfo, TRUE);

	/* Start decompressor */
	jpeg_start_decompress(&cinfo);

	printf("output_width=%d, output_height=%d, output_components=%d\n",
		cinfo.output_width,
		cinfo.output_height,
		cinfo.output_components);

	image->x=cinfo.output_width;
	image->y=cinfo.output_height;
	image->depth=cinfo.output_components;

	scanline_len = cinfo.output_width * cinfo.output_components;
	image->pixels=malloc(cinfo.output_width * cinfo.output_height * cinfo.output_components);

	while (scanline_count < cinfo.output_height) {
		output_data = (image->pixels + (scanline_count * scanline_len));
		jpeg_read_scanlines(&cinfo, &output_data, 1);
		scanline_count++;
	}

	/* Finish decompressing */
	jpeg_finish_decompress(&cinfo);

	jpeg_destroy_decompress(&cinfo);

	fclose(fff);

	return 0;
}

static int store_jpeg(char *filename, struct image_t *image) {

	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	int quality=90; /* % */
	int i;

	FILE *fff;

	JSAMPROW row_pointer[1];
	int row_stride;

	/* setup error handler */
	cinfo.err = jpeg_std_error(&jerr);

	/* initialize jpeg compression object */
	jpeg_create_compress(&cinfo);

	/* Open file */
	fff = fopen(filename, "wb");
	if (fff==NULL) {
		fprintf(stderr, "can't open %s: %s\n",
			filename,strerror(errno));
		return -1;
	}

	jpeg_stdio_dest(&cinfo, fff);

	/* Set compression parameters */
	cinfo.image_width = image->x;
	cinfo.image_height = image->y;
	cinfo.input_components = image->depth;
	cinfo.in_color_space = JCS_RGB;
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);

	/* start compressing */
	jpeg_start_compress(&cinfo, TRUE);

	row_stride=image->x*image->depth;

	for(i=0;i<image->y;i++) {
		row_pointer[0] = & image->pixels[i * row_stride];
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	/* finish compressing */
	jpeg_finish_compress(&cinfo);

	/* close file */
	fclose(fff);

	/* clean up */
	jpeg_destroy_compress(&cinfo);

	return 0;
}

int main(int argc, char **argv) {

	struct image_t image;

	/* Check command line usage */
	if (argc<2) {
		fprintf(stderr,"Usage: %s image_file\n",argv[0]);
		return -1;
	}

	/* Load an image */
	load_jpeg(argv[1],&image);

	Display *display = XOpenDisplay(NULL);
	int screen_num = DefaultScreen(display);
	Window root = RootWindow(display,screen_num);
	Visual *visual = DefaultVisual(display,screen_num);

	//char *data = (char*)malloc(256*256*4);

	XImage *img = XCreateImage(display,visual,image.depth,XYPixmap,
		0,image.pixels,256,256,32,0);
	fprintf(stdout, "Sucessfully create image\n");
	/*	
	int count = 256 * 256;
	for (int i = 0; i < count; ++i)
	{
		fprintf(stdout, "%d\n", i);
		XAddPixel(img,(long)image.pixels[i]);
	}
	fprintf(stdout,"end\n");
	*/

	Window win = XCreateSimpleWindow(display,root,50,50,256,256,1,0,0);
	XSelectInput(display,win,ExposureMask|KeyPressMask);
	XMapWindow(display,win);
	XPutImage(display,win,DefaultGC(display,screen_num),img,0,0,0,0,256,256);

	XEvent event;
	while(1)
	{
		XNextEvent(display,&event);
	}
	return 0;
	/* Write data back out to disk */
    //store_jpeg("out.jpg",&image);

	return 0;
}
