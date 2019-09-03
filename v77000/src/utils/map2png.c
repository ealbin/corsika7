#define png_infopp_NULL (png_infopp)NULL 
#define int_p_NULL (int*)NULL

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <png.h>

/**********************************************************
 ********************* M A P 2 P N G **********************
 **********************************************************

 Transforms particle track maps as written out by the air
 shower simulation program CORSIKA in the 'PLOTSH2' option.

 To compile the program, run (you need libpng.so):

 gcc [-L<path to libpng>] -lm -lpng map2png.c -o map2png


 This program was written by Fabian Schmidt, University of Leeds, UK
 (C) Fabian Schmidt, July 2005
 **********************************************************
 **********************************************************/


/* where map is copied and transformed to */
float* picture;
/* PNG arrays */
png_uint_16p pic16;  // if 16bit color depth
png_bytep pic8;      // if 8bit color depth 
unsigned int NpixX = 0, NpixY = 0;

/* settings */
char *maproot, *outputfile;
char projstr[] = "xz";
int bit_depth = 8;
int donb = 0;

/* track numbers -> color */
int logscale = 0;  // 0 -> linear scale
float emR = 1., emG = 0., emB = 0.;
float muR = 0., muG = 1., muB = 0.;
float hdR = 0., hdG = 0., hdB = 1.;

/* background color */
float bgR = 0., bgG = 0., bgB = 0.;


float findmaxvalue()
{
  unsigned int N = 3*NpixX * NpixY;
  float max = picture[0];
  unsigned int i = 0;
  for (i=1; i < N; i++) {
    if (picture[i] > max)  max = picture[i];
  }

  return max;
}

float findminvalue()
{
  unsigned int N = 3*NpixX * NpixY;
  float min = picture[0];
  unsigned int i = 0;
  for (i=1; i < N; i++) {
    if (picture[i] < min)  min = picture[i];
  }

  return min;
}

void reset_picture()
{
  unsigned int N = NpixX * NpixY;
  unsigned int i = 0;
  float *p = picture;
  for (i=0; i < N; i++) {
    *p++ = bgR;
    *p++ = bgG;
    *p++ = bgB;
  }
}

void scalepicture(float colormax)
{
  float fmin = findminvalue(), fmax = findmaxvalue(), factor=0.;
  if (fmin == fmax)  factor = 1.;
  else 
    factor = colormax / (fmax - fmin);

  printf("Maximal/minimal pixel value in map: %.4f/%.4f; rescaling picture by %.4f.\n",
	 fmax, fmin, factor);

  unsigned int N = 3*NpixX * NpixY;
  unsigned int i = 0;
  float *p = picture;
  for (i=0; i < N; i++) {
    *p = factor * (*p - fmin);  
    *p++;
  }
}


int read_map(FILE *fp, float R, float G, float B)
{
  unsigned int i,j;
  float *p = picture;
  float *b, *buffer = (float*) calloc(NpixX+2, sizeof(float));
  if (buffer == NULL) {
    printf("Couldn't allocate memory.\n");
    return -1;
  }

  short ir = R >= 0. ? 0 : 1,
    ig = G >= 0. ? 0 : 1,
    ib = B >= 0. ? 0 : 1;
  if (logscale) {
    if (ir) R = fabs(R);
    if (ig) G = fabs(G);
    if (ib) B = fabs(B);
  }
	
  for (j = 0; j < NpixY; j++) {
    if (fread(buffer, 4, NpixX+2, fp) != NpixX+2) {
      printf("Row %d: reading from file failed.\n", j);
      free(buffer);
      return -1;
    }
    b = buffer + 1;
/*     fill pixels */
    if (logscale) {
      for (i=0; i < NpixX; i++) {
	*p++ += ir ? -log10(R * (*b) + 1.) : log10(R * (*b) + 1.);
	*p++ += ig ? -log10(G * (*b) + 1.) : log10(G * (*b) + 1.);
	*p++ += ib ? -log10(B * (*b++) + 1.) : log10(B * (*b++) + 1.);
      }
    } 
    else {
      for (i=0; i < NpixX; i++) {
	*p++ += R * (*b); 
	*p++ += G * (*b); 
	*p++ += B * (*b++);
      }
    } 
  }

  free(buffer);
  return 0;
}

void write_png(const char *file_name, png_uint_32 width, png_uint_32 height,
               int bit_depth );


void print_usage()
{
  printf("\n Usage: map2png [options] <CORSIKA DATnnnnnn file> <output PNG file>\n\n");
  printf(" Options: (options separated by '/' are mutually exclusive)\n");
  printf(" -xy / -xz / -yz       Projection to be used (default: -xz).\n");  
  printf(" -lin / -log           Linear / log color scale (default: -lin).\n");
  printf(" -depth d              Image bit depth; currently only depth 8 supported.\n");
  printf(" -bg r,g,b             Background color (default: 0,0,0)\n");
  printf(" -em r,g,b             Color used for electromagnetic tracks (default: 1.,0,0)\n");
  printf(" -mu r,g,b             Color used for muonic tracks (default: 0,1,0)\n");
  printf(" -hd r,g,b             Color used for hadronic tracks (default: 0,0,1)\n");
  printf("\n Colors must be given as three float numbers separated by commas.\n");
  printf(" If you want dark tracks on bright background, you should specifiy\n");
  printf(" NEGATIVE colors for particle species (or invert the picture afterwards).\n");
  printf(" Image will be scaled to appropriate color range before writing.\n");
  printf(" map2png will read DATnnnnnn.<particle type>_<projection>.map\n for the three particle types.\n");

  printf("\n");
}

int parse_args(int argc, char **argv)
{
  int i, fncount = 0;

  for (i = 1; i < argc; i++) {
    if (*argv[i] == '-') {  // option
      if (!strncmp(argv[i]+1, "xy", 2) || 
	  !strncmp(argv[i]+1, "xz", 2) || 
	  !strncmp(argv[i]+1, "yz", 2) ) {
	strncpy(projstr, argv[i]+1, 2);
	printf("Using %s-projection.\n", projstr);
      }
      else if (!strncmp(argv[i]+1, "lin", 3)) {
	logscale = 0;
      }
      else if (!strncmp(argv[i]+1, "log", 3)) {
	logscale = 1;
      }
      else if (!strncmp(argv[i]+1, "depth", 5)) {
/* 	only 8 bit supported currently */
	if (atoi(argv[++i]) != 8) {
	  printf("Only 8 bit color depth currently supported.\n");
	  return 1;
	}
      }
      else if (!strncmp(argv[i]+1, "bg", 2)) {
	if (sscanf(argv[++i], " %f,%f,%f", &bgR, &bgG, &bgB) != 3) {
	  printf("Couldn't decipher color specification: %s\n", argv[i]);
	  return 1;
	}
	else 
	  printf("Background color: %f/%f/%f\n", bgR, bgG, bgB);
      }
      else if (!strncmp(argv[i]+1, "em", 2)) {
	if (sscanf(argv[++i], " %f,%f,%f", &emR, &emG, &emB) != 3) {
	  printf("Couldn't decipher color specification: %s\n", argv[i]);
	  return 1;
	}
	else 
	  printf("EM track color: %f/%f/%f\n", emR, emG, emB);
      }
      else if (!strncmp(argv[i]+1, "mu", 2)) {
	if (sscanf(argv[++i], " %f,%f,%f", &muR, &muG, &muB) != 3) {
	  printf("Couldn't decipher color specification: %s\n", argv[i]);
	  return 1;
	}
	else 
	  printf("Muon track color: %f/%f/%f\n", muR, muG, muB);
      }
      else if (!strncmp(argv[i]+1, "hd", 2)) {
	if (sscanf(argv[++i], " %f,%f,%f", &hdR, &hdG, &hdB) != 3) {
	  printf("Couldn't decipher color specification: %s\n", argv[i]);
	  return 1;
	}
	else 
	  printf("Hadronic track color: %f/%f/%f\n", hdR, hdG, hdB);
      }
      else {
	printf("Unknown option: %s\n", argv[i]);
	return 1;
      }
    }
    else {  /* filename */
      if (!fncount)  maproot = argv[i];
      else if (fncount == 1)  outputfile = argv[i];
      else 
	printf("Warning: unused argument '%s'.\n", argv[i]);
      fncount++;
    }
  }

  return 0;
}


int main(int argc, char **argv)
{
  if (argc < 3)  {
    print_usage();
    return 0;
  }

  if (parse_args(argc, argv))  return -1;

  if (sizeof(float) != 4) {
    printf("WARNING: size of float != 4 ! Adapt source code to your architecture.\n");
    return 0;
  }

  char mapfile[500];
  sprintf(mapfile, "%s.em_%s.map", maproot, projstr);
  FILE* fpem = fopen(mapfile, "r");
  if (!fpem) {
    printf("\n Couldn't open file '%s' for reading.\n", mapfile);
  }
  sprintf(mapfile, "%s.mu_%s.map", maproot, projstr);
  FILE* fpmu = fopen(mapfile, "r");
  if (!fpmu) {
    printf("\n Couldn't open file '%s' for reading.\n", mapfile);
  }
  sprintf(mapfile, "%s.hd_%s.map", maproot, projstr);
  FILE* fphd = fopen(mapfile, "r");
  if (!fphd) {
    printf("\n Couldn't open file '%s' for reading.\n", mapfile);
  }
  if (!fpem && !fpmu && !fphd)  return -1;

/*   read 'header' */
  float buf[4];
  if (fpem && fread(&buf[0], 4, 4, fpem) != 4)  {
    printf("Couldn't read from 'em' file.\n");
    return -1;
  }
  NpixX = (unsigned int) buf[1];
  NpixY = (unsigned int) buf[2];
  if (fpmu && fread(&buf[0], 4, 4, fpmu) != 4)  {
    printf("Couldn't read from 'mu' file.\n");
    return -1;
  }
  if ((unsigned int) buf[1] != NpixX || (unsigned int) buf[2] != NpixY) {
    printf("Map dimensions differ ! Cannot proceed.\n");
    return -1;
  }
  if (fphd && fread(&buf[0], 4, 4, fphd) != 4)  {
    printf("Couldn't read from 'hd' file.\n");
    return -1;
  }
  if ((unsigned int) buf[1] != NpixX || (unsigned int) buf[2] != NpixY) {
    printf("Map dimensions differ ! Cannot proceed.\n");
    return -1;
  }

/*   allocate picture & buffer */
  printf("Reading map of %d x %d pixels.\n", NpixX, NpixY);
  picture = (float*) calloc(3*NpixX*NpixY, sizeof(float));
  if (picture == NULL) {
    printf("Couldn't allocate memory for picture !\n");
    return -1;
  }
  reset_picture();

/*   Convert track map to 'float' color picture */
  if (fpem && read_map(fpem, emR, emG, emB)) {
    printf("Error while reading from em map file.\n");
    return -1;
  }
  if (fpmu && read_map(fpmu, muR, muG, muB)) {
    printf("Error while reading from mu map file.\n");
    return -1;
  }
  if (fphd && read_map(fphd, hdR, hdG, hdB)) {
    printf("Error while reading from hd map file.\n");
    return -1;
  }
  if (fpem)  fclose(fpem);    
  if (fpmu)  fclose(fpmu);    
  if (fphd)  fclose(fphd);    

/*   Convert float picture to PNG array */
  if (bit_depth == 8) {
     unsigned int i, N = 3*NpixX*NpixY;
     pic8 = (png_bytep) calloc(N, sizeof(png_byte));
     if (pic8 == NULL) {
       printf("Couldn't allocate memory for picture !\n");
       free(picture);
       return -1;
     }
/*      Scale picture to dynamic range available */
     scalepicture(255.);

     png_bytep p8 = pic8;
     float *pf = picture;
     for (i = 0; i < N; i++)
       *p8++ = (png_byte) *pf++;
  }
  else {

  }
  free(picture);

/*   Now write PNG */
  printf("Writing PNG picture to '%s'.\n", outputfile);
  write_png(outputfile, NpixX, NpixY, bit_depth);
  if (bit_depth == 8)  free(pic8);
  return 0;
}


/* write PNG file */
/* bit_depth is either 8 or 16 (function chooses picture pointer accordingly) 
*/
void write_png(const char *file_name, png_uint_32 width, png_uint_32 height,
               int bit_depth )
 
{
   FILE *fp;
   png_structp png_ptr;
   png_infop info_ptr;
/*    png_colorp palette; */

   /* open the file */
   fp = fopen(file_name, "wb");
   if (fp == NULL) {
     printf("Couldn't open PNG file '%s' for writing.\n", file_name);
     return;
   }

   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also check that
    * the library version is compatible with the one used at compile time,
    * in case we are using dynamically linked libraries.  REQUIRED.
    */
   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                     NULL, NULL, NULL);

   if (png_ptr == NULL)
   {
      fclose(fp);
      return;
   }

   /* Allocate/initialize the image information data.  REQUIRED */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL)
   {
      fclose(fp);
      png_destroy_write_struct(&png_ptr,  png_infopp_NULL);
      return;
   }

   /* Set error handling.  REQUIRED if you aren't supplying your own
    * error handling functions in the png_create_write_struct() call.
    */
   if (setjmp(png_jmpbuf(png_ptr)))
   {
      /* If we get here, we had a problem reading the file */
      fclose(fp);
      png_destroy_write_struct(&png_ptr, &info_ptr);
      return;
   }

   /* set up the output control if you are using standard C streams */
   png_init_io(png_ptr, fp);

   /* Set the image information here.  Width and height are up to 2^31,
    * bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
    * the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
    * PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
    * or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
    * PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
    * currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
    */
   png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, 
                PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, 
                PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

   /* optional significant bit chunk */
   /* otherwise, if we are dealing with a color image then */
   /*   sig_bit.red = true_red_bit_depth;
   sig_bit.green = true_green_bit_depth;
   sig_bit.blue = true_blue_bit_depth;
   png_set_sBIT(png_ptr, info_ptr, sig_bit);*/


   /* Optionally write comments into the image */
   /*   text_ptr[0].key = "Title";
   text_ptr[0].text = "Mona Lisa";
   text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
   text_ptr[1].key = "Author";
   text_ptr[1].text = "Leonardo DaVinci";
   text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
   text_ptr[2].key = "Description";
   text_ptr[2].text = "<long text>";
   text_ptr[2].compression = PNG_TEXT_COMPRESSION_zTXt;
   png_set_text(png_ptr, info_ptr, text_ptr, 3);
   */

   /* Write the file header information.  REQUIRED */
   png_write_info(png_ptr, info_ptr);

   /* The easiest way to write the image (you may have a different memory
    * layout, however, so choose what fits your needs best).  You need to
    * use the first method if you aren't handling interlacing yourself.
    */
   png_uint_32 k;
   png_bytep row_pointers[height];

   if (bit_depth == 8) {
     for (k = 0; k < height; k++)
       row_pointers[k] = pic8 + k * width * 3;
   }
   else if (bit_depth == 16) {
     for (k = 0; k < height; k++)
       row_pointers[k] = (png_bytep) (pic16 + k * width * 3);
   }
/* write out the entire image data in one call */
   png_write_image(png_ptr, row_pointers);

   /* It is REQUIRED to call this to finish writing the rest of the file */
   png_write_end(png_ptr, info_ptr);

   /* clean up after the write, and free any memory allocated */
   png_destroy_write_struct(&png_ptr, &info_ptr);

   /* close the file */
   fclose(fp);

   return;
}


