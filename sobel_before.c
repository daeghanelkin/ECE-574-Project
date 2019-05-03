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

/* Filters */
static int sobel_x_filter[3][3]={{-1,0,+1},{-2,0,+2},{-1,0,+1}};
static int sobel_y_filter[3][3]={{-1,-2,-1},{0,0,0},{1,2,+1}};

/* Structure describing the image */
struct image_t {
    int x;
    int y;
    int depth;    /* bytes */
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
    int val;
    
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
    
    for(y=ystart;y<yend;y++) {
        for(x=1;x<old->x-1;x++) {
            val=0;
            for(d=0; d<3; d++) {
                sum=0;
                for(k=-1;k<2;k++) {
                    for(l=-1;l<2;l++) {
                    color=old->pixels[((y+l)*width)+(x*depth+d+k*depth)];
                    sum+=color * (*filter)[k+1][l+1];
                    }
                }

                if (sum<0) sum=0;
                if (sum>255) sum=255;
                
                new->pixels[(y*width)+x*depth+d]=sum;
                val |= (sum<<8*(2-d));
            }
            // Update the pixel at x,y with the new value we have
            XPutPixel(image->img, x, y, (long)val);
        }
        XPutImage(image->dpy,image->win,DefaultGC(image->dpy,image->screen),
                   image->img,0,0,0,0,old->x,old->y);
    } 
    
    return NULL;
}

// Combine routine
static int combine(struct image_t *s_x,
            struct image_t *s_y,
            struct image_t *new,
            struct XImage_data_t *image) {
    int out;
    int val;
    int pixel;
    int x,y,d;
    int xsize = s_x->x;
    int ysize = s_x->y;
    int depth = s_x->depth;

    for(y=0;y<ysize;y++) {
        for(x=0; x < xsize; x++){
            val = 0;
            for(d=0; d<3; d++) {
                pixel = (y*xsize*depth)+x*depth+d;
                out=sqrt((s_x->pixels[pixel]*s_x->pixels[pixel])
                        +(s_y->pixels[pixel]*s_y->pixels[pixel]));
                if (out>255) out=255;
                if (out<0) out=0;
                new->pixels[pixel]=out;
                val |= (out<<8*(2-d));
            }

            // Update the pixel at x,y with the new value we have
            XPutPixel(image->img, x, y, (long)val);
        }
        XPutImage(image->dpy,image->win,DefaultGC(image->dpy,image->screen),
                  image->img,0,0,0,0,s_x->x,s_x->y);
    }

    return 0;
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

int main(int argc, char **argv) {

    struct image_t image, sobel_x, sobel_y, new_image;
    struct convolve_data_t sobel_data;
    struct XImage_data_t x_data;
    int val, x, y, d;
    int xsize, ysize, depth;
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

    /* Load an image */
    load_jpeg(argv[1],&image);

    /* Allocate space for output image */
    new_image.x=image.x;
    new_image.y=image.y;
    new_image.depth=image.depth;
    new_image.pixels=malloc(image.x*image.y*image.depth*sizeof(char));

    /* Allocate space for output image */
    sobel_x.x=image.x;
    sobel_x.y=image.y;
    sobel_x.depth=image.depth;
    sobel_x.pixels=malloc(image.x*image.y*image.depth*sizeof(char));

    /* Allocate space for output image */
    sobel_y.x=image.x;
    sobel_y.y=image.y;
    sobel_y.depth=image.depth;
    sobel_y.pixels=malloc(image.x*image.y*image.depth*sizeof(char));

    xsize = image.x;
    ysize = image.y;
    depth = image.depth;

    // Initial X11 window setup
    display = XOpenDisplay(NULL);
    screen_num = DefaultScreen(display);
    root = RootWindow(display,screen_num);
    visual = DefaultVisual(display,screen_num);
    
    // Create X11 Image using input image parameters
    // Although we load images with 24bpp, X11 only has 8, 16, and 32bpp 
    // capabilities so we use 32bpp
    data = (char *)malloc(xsize*ysize*4);
    img = XCreateImage(display,visual,DefaultDepth(display,screen_num),
                       ZPixmap,0,data,xsize,ysize,32,0);

    // Setup window for displaying
    win = XCreateSimpleWindow(display,root,50,50,xsize,ysize,1,0,0);
    XMapWindow(display,win);

    // Setup X11 image data struct
    x_data.dpy = display;
    x_data.img = img;
    x_data.screen = screen_num;
    x_data.win = win;

    // Display original image to X11 Screen
    for(y = 0; y < ysize; y++) {
        for(x = 0; x < xsize; x++) {
            val = 0;
            for(d = 0; d<3; d++) {
                val |= (image.pixels[(y*xsize*depth)+x*depth+d]<<8*(2-d));
            }
            XPutPixel(img, x, y, (long)val);
        }
    }
    // Set window name so you know what's running
    XStoreName(display, win, "Sobel X Convolution");
    sobel_data.old=&image;
    sobel_data.new=&sobel_x;
    sobel_data.filter=&sobel_x_filter;
    sobel_data.ystart=0;
    sobel_data.yend=image.y;
    generic_convolve((void *)&sobel_data, &x_data);

    // Set window name so you know what's running
    XStoreName(display, win, "Sobel Y Convolution");
    sobel_data.old=&image;
    sobel_data.new=&sobel_y;
    sobel_data.filter=&sobel_y_filter;
    sobel_data.ystart=0;
    sobel_data.yend=image.y;
    generic_convolve((void *)&sobel_data, &x_data);

    // Set window name so you know what's running
    XStoreName(display, win, "Combine");
    combine(&sobel_x, &sobel_y, &new_image, &x_data);
    // Destroy image now that we're done
    XDestroyImage(img);

    // Store image
    store_jpeg("out.jpg",&new_image);

    return 0;
}
