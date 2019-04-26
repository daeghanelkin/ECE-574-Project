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

// Structure for convolve data
struct convolve_data_t {
	struct image_t *old;
	struct image_t *new;
	int (*filter)[3][3];
	int ystart;
	int yend;
};

// Structure for X11 window information
struct XImage_data_t {
    XImage *img;
    Display *dpy;
    int screen;
    Window win;
};

/* very inefficient convolve code */
static void *generic_convolve(void *argument) {

	int x,y,k,l,d;
	uint32_t color;
	int sum,depth,width;

	struct image_t *old;
	struct image_t *new;
	int (*filter)[3][3];
	struct convolve_data_t *data;
	int ystart, yend;

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

    for(x=1;x<old->x-1;x++) {
        for(y=ystart;y<yend;y++) {
            for(d=0;d<3;d++) {
                sum=0;
				for(k=-1;k<2;k++) {
                    for(l=-1;l<2;l++) {
                        color=old->pixels[((y+l)*width)+(x*depth+d+k*depth)];
						sum+=color*(*filter)[k+1][l+1];
                    }
				}

				if (sum<0) sum=0;
				if (sum>255) sum=255;

				new->pixels[((y-data->ystart)*width)+x*depth+d]=sum;
            }
        }
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

static int combine(struct image_t *s_x, struct image_t *s_y, struct image_t *new, struct XImage_data_t *image) {
	int out;
    int val;
    int x,y,d;
    int xsize = s_x->x;
    int ysize = s_x->y;
    int depth = s_x->depth;

    for(y=0;y<ysize;y++) {
        for(x=0; x < xsize; x++){
            val = 0;
            for(d=0; d<3; d++) {
                out=sqrt((s_x->pixels[(y*xsize*depth)+x*depth+d]*s_x->pixels[(y*xsize*depth)+x*depth+d])
                        +(s_y->pixels[(y*xsize*depth)+x*depth+d]*s_y->pixels[(y*xsize*depth)+x*depth+d]));

                if (out>255) out=255;
                if (out<0) out=0;

                new->pixels[(y*xsize*depth)+x*depth+d]=out;
                val |= (out<<8*(2-d));
            }
            XPutPixel(image->img, x, y, (long)val);
        }
        XPutImage(image->dpy,image->win,DefaultGC(image->dpy,image->screen),image->img,0,0,0,0,s_x->x,s_x->y);
	}

	return 0;
}

int main(int argc, char **argv) {

	struct image_t image, sobel_x, sobel_y, new_image, tail;
	struct convolve_data_t sobel_data;
    struct XImage_data_t x_data;
	double start_time=0, load_time=0, store_time=0, convolve_time=0, combine_time=0;
	int result, rank, numtasks,y,x,d;
	int image_info[3];		//x,y,depth
	int i, j;
    int tailstart, tailsize;
	long val = 0;
	char *data;
    Display *display;
    int screen_num;
    Window root, win;
    Visual *visual;
    XImage *img;
    XEvent event;
    int *counts;
    int *displacements;

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
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
    /* and this processes' rank is */
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // Initial setup on base node
	if (rank == 0) {
		printf("Size: %d\nRank: %d\n",numtasks, rank);

		start_time=MPI_Wtime();
		load_jpeg(argv[1],&image);		//load the image
		load_time=MPI_Wtime();				
		printf("Load time: %lf\n",load_time-start_time);

		image_info[0] = image.x;			//horizontal pixel #
		image_info[1] = image.y;			//vertical pixel #
		image_info[2] = image.depth;	//image depth

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
        XSelectInput(display, win, ExposureMask | KeyPressMask);
        XStoreName(display, win, "Base Image");

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
		XPutImage(display,win,DefaultGC(display,screen_num),img,0,0,0,0,image.x,image.y);

        // Allocate space for rank counts (needed for MPI_Gatherv)
        counts=calloc(numtasks, sizeof(int));
        for(i=0; i<numtasks; i++){
            counts[i]=image.x*image.depth;
        }

        // Allocate space for displacements (needed for MPI_Gatherv)
        displacements=calloc(numtasks, sizeof(int));
	}

    // Make all pcoesses wait
	MPI_Barrier(MPI_COMM_WORLD);		

    // Broadcast paremeters
	result = MPI_Bcast(image_info, 3,MPI_INT, 0, MPI_COMM_WORLD);
    if (result != MPI_SUCCESS) {
        printf("Error broadcasting image info.\n");
        MPI_Abort(MPI_COMM_WORLD, result);
    }

    // Set up image for remaining ranks
	if (rank != 0) {
		image.x = image_info[0];			//set image width
		image.y = image_info[1];			//set image heught
		image.depth = image_info[2];	//set image depth
		image.pixels=calloc(image.x*image.y*image.depth,sizeof(char));	//allocate memory
	}

	// Allocate space for final image
	new_image.x = image.x;
	new_image.y = image.y;
	new_image.depth = image.depth;
	new_image.pixels=calloc(image.x*image.y*image.depth,sizeof(char));

	// Allocate space for Sobel X image
	sobel_x.x = image.x;
	sobel_x.y = image.y;
	sobel_x.depth = image.depth;
	sobel_x.pixels=calloc(image.x*image.y*image.depth,sizeof(char));

	// Allocate space for Sobel Y image
	sobel_y.x = image.x;
	sobel_y.y = image.y;
	sobel_y.depth = image.depth;
	sobel_y.pixels=calloc(image.x*image.y*image.depth,sizeof(char));

    // Allocate space for tail image
    if (rank == 0) {
        tail.x = image.x;
        tail.y = image.y;
        tail.depth = image.depth;
        tail.pixels = calloc(image.x*image.y*image.depth,sizeof(char));
    }
		
	// Broadcast the input pixels for sobel convoluctions
	result = MPI_Bcast(image.pixels, image.x*image.y*image.depth, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (result != MPI_SUCCESS) {
        printf("Error broadcasting image.\n");
        MPI_Abort(MPI_COMM_WORLD, result);
    }

    // Setup data for sobel_x convolution
    sobel_data.old = &image;
    sobel_data.new = &new_image;
    sobel_data.filter = &sobel_x_filter;
    sobel_data.ystart = image.y / numtasks * rank;
    sobel_data.yend = image.y / numtasks * rank + 1;

    if (rank == 0) {
        XStoreName(display, win, "Sobel X Filter");
    }

    // Loop through the amount of even lines available
    for(i = 0; i < image.y / numtasks; i++) {
        generic_convolve((void *)&sobel_data);

        // Adjust displacements for putting lines into sobel buffer
        if (rank == 0) {
            for(j = 0; j < numtasks; j++){
                displacements[j]=(j*image.y/numtasks*image.x*image.depth)+(i*image.x*image.depth);
            }
        }

        // Gather a line at a time from each rank
        MPI_Gatherv(new_image.pixels,
                   image.x*image.depth,
                   MPI_CHAR,
                   sobel_x.pixels,
                   counts,
                   displacements,
                   MPI_CHAR,
                   0,
                   MPI_COMM_WORLD);

        
        // Update image
        if (rank == 0) {
            for(y = sobel_data.ystart; y < image.y; y+=image.y/numtasks) {
                for(x = 0; x < image.x; x++) {
                    val = 0;
                    for(d = 0; d<3; d++) {
                        val |= (sobel_x.pixels[(y*image.x*image.depth)+x*image.depth+d]<<8*(2-d));
                    }
                    XPutPixel(img, x, y, (long)val);
                }

                // Display image
                XPutImage(display,win,DefaultGC(display,screen_num),img,0,0,0,0,image.x,image.y);
            }
        }

        // Increase the start and stop y values for each rank
        if (sobel_data.yend < image.y/numtasks*(rank+1)) {
            sobel_data.ystart++;
            sobel_data.yend++;
        }
    }

    // Calculate tail end of image
    if (rank == 0 && image.y%numtasks) {
        sobel_data.old = &image;
        sobel_data.new = &tail;
        sobel_data.filter = &sobel_x_filter;
        sobel_data.ystart = (image.y/numtasks)*numtasks;
        sobel_data.yend = image.y;

        generic_convolve((void *)&sobel_data);

        tailstart = (image.y/numtasks)*numtasks;
        tailsize = image.y-(image.y/numtasks)*numtasks;
        tailsize *= image.x*image.depth;

        memcpy(&sobel_x.pixels[tailstart*image.x*image.depth], &tail.pixels[0], tailsize);
    }

    // Setup data for sobel_y convolution
    sobel_data.old = &image;
    sobel_data.new = &new_image;
    sobel_data.filter = &sobel_y_filter;
    sobel_data.ystart = image.y / numtasks * rank;
    sobel_data.yend = image.y / numtasks * rank + 1;

    if (rank == 0) {
        XStoreName(display, win, "Sobel Y Filter");
    }

    // Loop through the amount of even lines available
    for(i = 0; i < image.y / numtasks; i++) {
        generic_convolve((void *)&sobel_data);

        // Adjust displacements for putting lines into sobel buffer
        if (rank == 0) {
            for(j = 0; j < numtasks; j++){
                displacements[j]=(j*image.y/numtasks*image.x*image.depth)+(i*image.x*image.depth);
            }
        }

        // Gather a line at a time from each rank
        MPI_Gatherv(new_image.pixels,
                   image.x*image.depth,
                   MPI_CHAR,
                   sobel_y.pixels,
                   counts,
                   displacements,
                   MPI_CHAR,
                   0,
                   MPI_COMM_WORLD);

        
        // Update image
        if (rank == 0) {
            for(y = sobel_data.ystart; y < image.y; y+=image.y/numtasks) {
                for(x = 0; x < image.x; x++) {
                    val = 0;
                    for(d = 0; d<3; d++) {
                        val |= (sobel_y.pixels[(y*image.x*image.depth)+x*image.depth+d]<<8*(2-d));
                    }
                    XPutPixel(img, x, y, (long)val);
                }

                // Display image
                XPutImage(display,win,DefaultGC(display,screen_num),img,0,0,0,0,image.x,image.y);
            }
        }

        // Increase the start and stop y values for each rank
        if (sobel_data.yend < image.y/numtasks*(rank+1)) {
            sobel_data.ystart++;
            sobel_data.yend++;
        }
    }

    // Calculate tail end of image
    if (rank == 0 && image.y%numtasks) {
        sobel_data.old = &image;
        sobel_data.new = &tail;
        sobel_data.filter = &sobel_y_filter;
        sobel_data.ystart = (image.y/numtasks)*numtasks;
        sobel_data.yend = image.y;

        generic_convolve((void *)&sobel_data);

        tailstart = (image.y/numtasks)*numtasks;
        tailsize = image.y-(image.y/numtasks)*numtasks;
        tailsize *= image.x*image.depth;

        memcpy(&sobel_y.pixels[tailstart*image.x*image.depth], &tail.pixels[0], tailsize);
    }

	if(rank == 0) {
        convolve_time = MPI_Wtime();

        x_data.dpy = display;
        x_data.img = img;
        x_data.screen = screen_num;
        x_data.win = win;

		/* Combine to form output */
		combine(&sobel_x, &sobel_y, &new_image, &x_data);
		combine_time = MPI_Wtime();

		store_jpeg("out.jpg",&new_image);
		store_time=MPI_Wtime();

 		printf("Convolve time: %lf\n",convolve_time-load_time);
 		printf("Combine time: %lf\n",combine_time-convolve_time);
 		printf("Store time: %lf\n",store_time-combine_time);
		printf("Total time = %lf\n",store_time-start_time);
        
        // Display image until key pressed
        while(1) {
            XNextEvent(display, &event);

            if(event.type == KeyPress)
                break;
        }

        // Clean up X windows
        XDestroyWindow(display, win);
        XCloseDisplay(display);
	}

	MPI_Finalize();


	return 0;
}
