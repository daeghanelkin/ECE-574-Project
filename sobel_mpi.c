/* Example sobel code for ECE574 -- Fall 2015 */
/* By Vince Weaver <vincent.weaver@maine.edu> */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <math.h>

#include <jpeglib.h>

#include <mpi.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

/* Filters */
static int sobel_x_filter[3][3]={{-1,0,+1},{-2,0,+2},{-1,0,+1}};
static int sobel_y_filter[3][3]={{-1,-2,-1},{0,0,0},{1,2,+1}};

/* Structure describing the image */
struct image_t {
	int x;
	int y;
	int depth;	/* bytes */
	unsigned char *pixels;
};

struct convolve_data_t {
	struct image_t *old;
	struct image_t *new;
	int (*filter)[3][3];
	int ystart;
	int yend;
};

// Struct for X11 Image data
struct XImage_data_t {
        XImage *img;
        Display *dpy;
        int screen;
        Window win;
};

/* very inefficient convolve code */
static void *generic_convolve(void *argument, struct XImage_data_t *image) {

	int x,y,k,l,d;
	uint32_t color;
	int sum,depth,width;

	struct image_t *old;
	struct image_t *new;
	int (*filter)[3][3];
	struct convolve_data_t *data;
	int ystart, yend;
	//long val = 0;
	/* Convert from void pointer to the actual data type */
	data=(struct convolve_data_t *)argument;
	old=data->old;
	new=data->new;
	filter=data->filter;

	ystart=data->ystart;
	yend=data->yend;

	depth=old->depth;
	width=old->x*old->depth;

	if (ystart==0) ystart=1;
	if (yend==old->y) yend=old->y-1;
	if (yend==ystart) yend++;

	for(d=0;d<3;d++) {
	   for(x=1;x<old->x-1;x++) {
	     for(y=ystart;y<yend;y++) {
				sum=0;
				for(k=-1;k<2;k++) {
		   		for(l=-1;l<2;l++) {
						color=old->pixels[((y+l)*width)+(x*depth+d+k*depth)];
						sum+=color * (*filter)[k+1][l+1];
		   		}
				}

				if (sum<0) sum=0;
				if (sum>255) sum=255;
				//added offset to put at begining of the buffer
				new->pixels[((y-data->ystart)*width)+x*depth+d]=sum;
	//			val |= (sum<<8*(2-d));
	    }
                        // Update the pixel at x,y with the new value we have
          //              XPutPixel(image->img, x, y, (long)val);
	  }
//	XPutImage(image->dpy,image->win,DefaultGC(image->dpy,image->screen),image->img,0,0,0,0,old->x,old->y);
	}
	return NULL;
}


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

static int combine(struct image_t *s_x,
			struct image_t *s_y,
			struct image_t *new,
			struct XImage_data_t *image) {
	int out;
       // int val;
        int x,y,d;
        int xsize = s_x->x;
        int ysize = s_x->y;
        int depth = s_x->depth;

    for(y=0;y<ysize;y++) {
                for(x=0; x < xsize; x++){
                        //val = 0;
                        for(d=0; d<3; d++) {
                                out=sqrt((s_x->pixels[(y*xsize*depth)+x*depth+d]*s_x->pixels[(y*xsize*depth)+x*depth+d])
                                        +(s_y->pixels[(y*xsize*depth)+x*depth+d]*s_y->pixels[(y*xsize*depth)+x*depth+d]));
                                if (out>255) out=255;
                                if (out<0) out=0;
                                new->pixels[(y*xsize*depth)+x*depth+d]=out;
//                                val |= (out<<8*(2-d));
                        }

                        // Update the pixel at x,y with the new value we have
//                        XPutPixel(image->img, x, y, (long)val);
                }
//        XPutImage(image->dpy,image->win,DefaultGC(image->dpy,image->screen),image->img,0,0,0,0,s_x->x,s_x->y);
	}

	return 0;
}

int main(int argc, char **argv) {

	struct image_t image,sobel_x,sobel_y,new_image,another,sobel_xa;
	struct convolve_data_t sobel_data[2];
	struct XImage_data_t x_data;
	double start_time=0,load_time=0,store_time,convolve_time,combine_time;
	int result, myid, numprocs,y,x,d;
	int buff[3];		//x,y,depth
	int i;
	int xsize, ysize, depth;
	long val = 0;
	char *data;
        Display *display;
        int screen_num;
        Window root, win;
        Visual *visual;
        XImage *img;

	/* Check command line usage */
	if (argc<2) {
		fprintf(stderr,"Usage: %s image_file\n",argv[0]);
		return -1;
	}

	/* Initialize MPI */
	result = MPI_Init(&argc,&argv);
	if (result != MPI_SUCCESS) {
		printf ("Error starting MPI program!.\n");
		MPI_Abort(MPI_COMM_WORLD, result);
	}


  /* find out how big the SPMD world is */
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  /* and this processes' rank is */
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	if (myid == 0) {
		printf("Size: %d\nRank: %d\n",numprocs, myid);

		start_time=MPI_Wtime();
		load_jpeg(argv[1],&image);		//load the image
		load_time=MPI_Wtime();				
		printf("Load time: %lf\n",load_time-start_time);

		buff[0] = image.x;			//horizontal pixel #
		buff[1] = image.y;			//vertical pixel #
		buff[2] = image.depth;	//image depth

		xsize = image.x;
		ysize = image.y;
		depth = image.depth;

		// Initial X11 window setup
		display = XOpenDisplay(NULL);
		screen_num = DefaultScreen(display);
		root = RootWindow(display,screen_num);
		visual = DefaultVisual(display,screen_num);

		// Create X11 Image using input image parameters
		// Although we load images with 24bpp, X11 only has 8, 16, and 32bpp capabilities so we use 32bpp
		data = (char *)malloc(image.x*image.y*4);
		img = XCreateImage(display,visual,DefaultDepth(display,screen_num),ZPixmap,0,data,image.x,image.y,32,0);

		// Setup window for displaying
		win = XCreateSimpleWindow(display,root,50,50,image.x,image.y,1,0,0);
		XMapWindow(display,win);

		// Setup X11 image data struct
		x_data.dpy = display;
		x_data.img = img;
		x_data.screen = screen_num;
		x_data.win = win;

		// Display original image to X11 Screen
		for(y = 0; y < image.y; y++) {
			for(x = 0; x < image.x; x++) {
				val = 0;
				for(d = 0; d<3; d++) {
					val |= (image.pixels[(y*image.x*image.depth)+x*image.depth+d]<<8*(2-d));
				}
				XPutPixel(img, x, y, (long)val);
			}
		}
		XPutImage(x_data.dpy,x_data.win,DefaultGC(x_data.dpy,x_data.screen),x_data.img,0,0,0,0,image.x,image.y);
	}

	MPI_Barrier(MPI_COMM_WORLD);		//make all pcoesses wait
	MPI_Bcast(buff, 3,MPI_INT, 0, MPI_COMM_WORLD);		//broadcast paremeters
	if (myid != 0) {
		image.x = buff[0];			//set image width
		image.y = buff[1];			//set image heught
		image.depth = buff[2];	//set image depth
		image.pixels=calloc(image.x*image.y*image.depth,sizeof(char));	//allocate memory

	}
	//allocate memory for sobels and new 
	sobel_x.pixels=calloc(buff[0]*buff[1]*buff[2],sizeof(char));
	sobel_xa.pixels=calloc(buff[0]*buff[1]*buff[2],sizeof(char));
	sobel_y.pixels=calloc(buff[0]*buff[1]*buff[2],sizeof(char));
	another.pixels=calloc(buff[0]*buff[1]*buff[2],sizeof(char));

	//set new image parameters
	new_image.x = buff[0];
	new_image.y = buff[1];
	new_image.depth = buff[2];
	new_image.pixels=calloc(buff[0]*buff[1]*buff[2],sizeof(char));
		
	//broadcast the input pixels for sobel convoluctions
	MPI_Bcast(image.pixels, buff[0]*buff[1]*buff[2],MPI_CHAR, 0, MPI_COMM_WORLD);

	//initilize sobel_x intermediate paremeters
	another.x = buff[0];
	another.y = buff[1];
	another.depth = buff[2];

	//initilize sobel_x paremeters
	sobel_x.x = buff[0];
	sobel_x.y = buff[1];
	sobel_x.depth = buff[2];

	//initialize sobel_y parameters
	sobel_y.x = buff[0];
	sobel_y.y = buff[1];
	sobel_y.depth = buff[2];

	//setup data for sobel_x convolution
	sobel_data[0].old=&image;
	sobel_data[0].new=&another;
 	sobel_data[0].filter=&sobel_x_filter;
	sobel_data[0].ystart = myid - numprocs;
	int total = buff[1] / numprocs;
	int j = 0, k = 0, h = 0,v;
	int xa = 0, ya = 0;
//	printf("Thread %d Entered\n", myid);
	for(v = 0; v < total; v++) {
		for(i = 0; i < numprocs; i++) {
			if(i == myid) {
				another.pixels=calloc(buff[0]*buff[1]*buff[2],sizeof(char));
				sobel_data[0].new=&another;
				sobel_data[0].ystart += numprocs;
				sobel_data[0].yend = sobel_data[0].ystart + 1;
		//		printf("Thread %d Start %d End %d\n", myid, sobel_data[0].ystart, sobel_data[0].yend);
				generic_convolve((void *)&sobel_data[0],&x_data);
				MPI_Barrier(MPI_COMM_WORLD);		//make all pcoesses wait

				MPI_Gather(another.pixels,              //source
					sobel_x.depth*sobel_x.x*(sobel_x.y/(numprocs)),    //count
					MPI_CHAR,  //type
					sobel_x.pixels,  //recieve buffer
					sobel_x.depth * sobel_x.x*(sobel_x.y/(numprocs)),    //count
					MPI_CHAR,    //type
					0,   //source
					MPI_COMM_WORLD);
			}
		}
		//printf("It is: %d on thread %d\n", i, myid);
		MPI_Barrier(MPI_COMM_WORLD);		//make all pcoesses wait
		if(myid == 0) {
			// Display original image to X11 Screen
			for(x = 0; x < xsize; x++) {
				for(y = sobel_data[0].ystart; y < (sobel_data[0].ystart + numprocs); y++) {
					val = 0;
					for(d = 0; d<3; d++) {
						val |= (sobel_x.pixels[(y*xsize*depth)+x*depth+d]<<8*(2-d));
						sobel_xa.pixels[(y*xsize*depth)+x*depth+d] = sobel_x.pixels[(y*xsize*depth)+x*depth+d];
					}
					XPutPixel(img, x, y, (long)val);
				}
			}
			XPutImage(x_data.dpy,x_data.win,DefaultGC(x_data.dpy,x_data.screen),x_data.img,0,0,0,0,image.x,image.y);
			printf("Y start %d, Y end %d\n",sobel_data[0].ystart, sobel_data[0].ystart + numprocs);
		}
		MPI_Barrier(MPI_COMM_WORLD);		//make all pcoesses wait
	}
/*
			for(k = (sobel_data[0].ystart); k < (sobel_data[0].ystart + numprocs); k++) {
				xa = 0;
				for(j=(k*image.x*image.depth)+xa*image.depth; j < ((k+1)*image.x*image.depth)+xa*image.depth; j) {
//					new->pixels[((y-data->ystart)*width)+x*depth+d]=sum;
					val = 0;
					for(h = 0; h <image.depth; h++) {
						val |= sobel_x.pixels[j++] <<8*(2-h);
					}
					XPutPixel(img, xa++, ya, (long)val);	
				}
				XPutImage(x_data.dpy,x_data.win,DefaultGC(x_data.dpy,x_data.screen),x_data.img,0,0,0,0,image.x,image.y);
				ya++;
			}
*/
	if(myid == 0) {
		for(i; i<image.x*image.y*image.depth; i++) {
			sobel_x.pixels[i] = sobel_xa.pixels[i];
		}
		store_jpeg("sobel_x.jpg",&sobel_x);
	}
	//if last process go to end of the y
	//if ((image.y % numprocs) && (myid == (numprocs - 1))) sobel_data[0].yend = image.y;
	//setup data for sobel_x convolution
/*
	sobel_data[1].old=&image;
	sobel_data[1].new=&new_image;
 	sobel_data[1].filter=&sobel_y_filter;
 	sobel_data[1].ystart = (image.y/numprocs) * myid;
 	sobel_data[1].yend = image.y/numprocs*(myid+1);
	//if last process go to end of the y
	if ((image.y % numprocs) && (myid == (numprocs - 1))) sobel_data[1].yend = image.y;
	//sobel)y convolution
 	generic_convolve((void *)&sobel_data[1],&x_data);

	convolve_time=MPI_Wtime();

	//gather all of the sobel_x processes and stoare into sobel_x
  MPI_Gather(another.pixels,              //source
           sobel_x.depth*sobel_x.x*(sobel_x.y/(numprocs)),    //count
           MPI_CHAR,  //type
           sobel_x.pixels,  //recieve buffer
           sobel_x.depth * sobel_x.x*(sobel_x.y/(numprocs)),    //count
           MPI_CHAR,    //type
           0,   //source
           MPI_COMM_WORLD);
	//gather all of the sobel_y process and store into sobel_y
  MPI_Gather(new_image.pixels,              //source
           sobel_x.depth*sobel_x.x*(sobel_x.y/(numprocs)),    //count
           MPI_CHAR,  //type
           sobel_y.pixels,  //recieve buffer
           sobel_x.depth * sobel_x.x*(sobel_x.y/(numprocs)),    //count
           MPI_CHAR,    //type
           0,   //source
           MPI_COMM_WORLD);
*/
	if(myid == 0) {

		/* Combine to form output */
//		combine(&sobel_x,&sobel_y,&new_image,&x_data);
//		combine_time=MPI_Wtime();

//		store_jpeg("out.jpg",&new_image);
		store_time=MPI_Wtime();

 		printf("Convolve time: %lf\n",convolve_time-load_time);
 		printf("Combine time: %lf\n",combine_time-convolve_time);
 		printf("Store time: %lf\n",store_time-combine_time);
		printf("Total time = %lf\n",store_time-start_time);
	}

	MPI_Finalize();

	return 0;
}
