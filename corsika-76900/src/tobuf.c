//************************************************************************
// Original code by Juergen Oehlschlaeger to write DAT file using C code
// 14/11/2012
// Bug fixes by Konrad Bernloehr
// 26/02/2014
// extensions for ICECUBE FiFo by J. van Santen, U Wisconsin-Madison
// Febr. 2015 (by D. Heck)
//************************************************************************

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if HAVE_CONFIG_H
#include "config.h"
#endif
#if __ICECUBE2__
#include <sys/stat.h>
#include <errno.h>
#endif

FILE *fmpatap;
#if __ICECUBE2__
int is_pipe;
#endif
#if __CERENKOV__
FILE **fmcetap;
static int count_fmcetap = 0;
#endif

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
#if __THIN__ || __ICECUBE1__
 void fwritempatap__( int *maxbuf, int *nsubbl, float outvect[6552] ) {
#else
 void fwritempatap__( int *maxbuf, int *nsubbl, float outvect[5733] ) {
#endif
#else
#if __THIN__ || __ICECUBE1__
 void fwritempatap_( int *maxbuf, int *nsubbl, float outvect[6552] ) {
#else
 void fwritempatap_( int *maxbuf, int *nsubbl, float outvect[5733] ) {
#endif
#endif

  //  printf("input %d, %d",*maxbuf,*nsubbl);

   union fltint {
         float funi;
         int iuni;
   } intflt4;

   intflt4.iuni = *maxbuf * *nsubbl * 4; // record length 22932 wihtout thin or 26208 with thin.

#if !__ICECUBE1__
   //   printf("outvect: %d, %f",intflt4.iuni,outvect[0]);
#endif

   if ( fmpatap != NULL ) /* Only try to write if we have an output file */
   {
      fwrite(&intflt4.funi, sizeof(float), 1, fmpatap);
      fwrite(outvect, intflt4.iuni, 1, fmpatap);
      fwrite(&intflt4.funi, sizeof(float), 1, fmpatap);
      /* Does it make sense to continue if there is an error (e.g. disk full) ? */
      /* if ( ferror(fmaptap) ) ... */
   }

}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
#if __ICECUBE2__
  void fopenmpatap__( const char *corsika_name, int *ibl, int *compress, int *pipe ) {
#else
 void fopenmpatap__( const char *corsika_name, int *ibl ) {
#endif
#else
#if __ICECUBE2__
 void fopenmpatap_( const char *corsika_name, int *ibl, int *compress, int *pipe ) {
#else
 void fopenmpatap_( const char *corsika_name, int *ibl ) {
#endif
#endif

   unsigned long length= (unsigned long) (*ibl);
   const char *s = strchr(corsika_name,' ');
   if ( s != NULL && (s-corsika_name) < length )
      length = (unsigned long) (s-corsika_name);
#if __ICECUBE2__
   char *file_name=malloc(length+4);
#else
   char *file_name=malloc(length+1);
#endif

   strncpy(file_name,corsika_name,length);
   file_name[length]='\0';

#if __ICECUBE2__
   if (*compress) {
       strncpy(file_name+length, ".gz", 4);
       length += 3;
   }

   if (*pipe) {
       int err = mkfifo(file_name, S_IRUSR | S_IWUSR);
       if (err != 0) {
           fprintf(stderr, "Could not create %s: %s", file_name, strerror(err));
           exit(err);
       }
   }

   if (*compress) {
       char *command=malloc(length+14);
       snprintf(command, length+14, "gzip -9 -c > %s", file_name);
       fmpatap = popen(command, "w");
       is_pipe = 1;
       if (!fmpatap) {
           fprintf(stderr, "Could not open pipe to %s: %s", command, strerror(errno));
           exit(errno);
       }
   } else {
       fmpatap = fopen( file_name, "w");
       is_pipe = 0;
   }
#endif

   if ( length != (unsigned long) (*ibl) )
   {
      fprintf(stderr,"\n\nFile name in fopenmpatap() function reported to have length %d\n", *ibl);
      fprintf(stderr,"but actually is only of length %lu: '%s'\n\n", length, file_name);
   }

   //     printf("file name: %s end\n",file_name) ;

#if __ICECUBE2__
    free(file_name);
#else
   /* If the null device was given, it is more efficient to not open anything */
   fmpatap = NULL;
   if ( strcmp(file_name,"/dev/null") != 0 )
   {
      fmpatap = fopen( file_name, "w");
      if ( fmpatap == NULL )
         perror(file_name);
   }
#endif
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
void fclosempatap__() {
#else
void fclosempatap_() {
#endif

#if __ICECUBE2__
    if (!is_pipe)
        fclose(fmpatap); // close test version of particle data file.
    else
        pclose(fmpatap);
#else
   if ( fmpatap != NULL )
      fclose(fmpatap); // close test version of particle data file.
   fmpatap = NULL;
#endif
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __CERENKOV__
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
#if __THIN__
 void fwritemcetap__( int *maxbf2, int *nsubbl, float outvect[6552], int *icerfi ) {
#else
 void fwritemcetap__( int *maxbf2, int *nsubbl, float outvect[5733], int *icerfi ) {
#endif
#else
#if __THIN__
   void fwritemcetap_( int *maxbf2, int *nsubbl, float outvect[6552], int *icerfi ) {
#else
   void fwritemcetap_( int *maxbf2, int *nsubbl, float outvect[5733], int *icerfi ) {
#endif
#endif

  //  printf("input %d, %d",*maxbf2,*nsubbl);
   const int id = (*icerfi) - 1; // fortran index -> c index

   union fltint {
         float funi;
         int iuni;
   } intflt4;

   intflt4.iuni = *maxbf2 * *nsubbl * 4; // record length 22932 wihtout thin or 26208 with thin.

   //   printf("outvect: %d, %f",intflt4.iuni,outvect[0]);

   if (id >= count_fmcetap || id<0) {
     fprintf(stderr, "Trying to write to CERENKOV output file %i which is not open ", id);
     exit(1);
   }
   
   if ( fmcetap[id] != NULL ) /* Only try to write if we have an output file */
   {
      fwrite(&intflt4.funi, sizeof(float), 1, fmcetap[id]);
      fwrite(outvect, intflt4.iuni, 1, fmcetap[id]);
      fwrite(&intflt4.funi, sizeof(float), 1, fmcetap[id]);
      /* Does it make sense to continue if there is an error (e.g. disk full) ? */
      /* if ( ferror(fmcetap) ) ... */
   }

}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
   void fopenmcetap__( const char *corsika_name, int *ibl, int *icerfi ) {
#else
   void fopenmcetap_( const char *corsika_name, int *ibl, int *icerfi ) {
#endif

   const int id = (*icerfi) - 1; // fortran index -> c index

   unsigned long length= *ibl;
   const char *s = strchr(corsika_name,' ');
   if ( s != NULL && (s-corsika_name) < length )
      length = (unsigned long) (s-corsika_name);
   char *file_name=malloc(length+1);

   strncpy(file_name,corsika_name,length);
   file_name[length]='\0';
   if ( length != (unsigned long) (*ibl) )
   {
      fprintf(stderr,"\n\nFile name in fopenmcetap() function reported to have length %d\n", *ibl);
      fprintf(stderr,"but actually is only of length %lu: '%s'\n\n", length, file_name);
   }

   //     printf("file name: %s end\n",file_name) ;

   if (id<0) {
     fprintf(stderr, "Cannot work with negative Cherenkov telescops ID %i\n", id);
     exit(1);
   }

   //   fprintf(stderr, "OPENING ID=%i file=%s \n", id, corsika_name);
   
   if (id>=count_fmcetap) {
     FILE** newFILE = malloc((id+1) * sizeof(FILE*));
     int i = 0;
     for (i = 0; i < count_fmcetap; ++i) {
       newFILE[i] = fmcetap[i];         
     }
     if (count_fmcetap>0) {
       free(fmcetap);
     }
     fmcetap  = newFILE;
     count_fmcetap = id + 1;
   }
   
   /* If the null device was given, it is more efficient to not open anything */
   fmcetap[id] = NULL;
   if ( strcmp(file_name,"/dev/null") != 0 )
   {
      fmcetap[id] = fopen( file_name, "w");
      if ( fmcetap[id] == NULL )
         perror(file_name);
   }
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
void fclosemcetap__(int *icerfi) {
#else
void fclosemcetap_(int *icerfi) {
#endif

  const int id = (*icerfi) - 1;
  
  if ((int)id >= count_fmcetap || id<0) {
     printf("Trying to close CERENKOV output file %d which is not open \n", id);
     exit(1);
   }

   if ( fmcetap[id] != NULL )
      fclose(fmcetap[id]); 
   fmcetap[id] = NULL;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#endif

#if __ICECUBE1__
// output buffer in C, where we have malloc()
// routines for ICECUBE1 FiFo by J. van Santen, U Wisconsin-Madison

typedef struct {
    size_t size, capacity;
    float *buffer;
} buffer_t;

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
void tobuf_c__(float *block, int *maxbuf, int *ifl, int *zeroes)
#else
void tobuf_c_(float *block, int *maxbuf, int *ifl, int *zeroes)
#endif
{
    static buffer_t writebuf = {0, 0, NULL};
    static int superblock = 21;
    size_t record_size = (*maxbuf)*sizeof(float);
    int i, j, imax;

    if (*ifl == 2) {
        // Just write the RUNE block
        memset(&block[(*maxbuf)], 0, (superblock-1)*(*maxbuf)*sizeof(float));
        fwritempatap_(maxbuf, &superblock, block);
    } else {
        // Expand the write buffer if needed (always in chunks of the
        // superblock size)
        if (writebuf.size == writebuf.capacity) {
            writebuf.capacity += superblock*(*maxbuf);
            writebuf.buffer = realloc(writebuf.buffer,
                writebuf.capacity*sizeof(float));
            if (writebuf.buffer == NULL) {
                fprintf(stderr, "Couldn't allocate %zu bytes for output buffer!\n",
                           writebuf.capacity*sizeof(float));
           //     exit(ENOMEM);
                  exit(-1);
            }
            memset(writebuf.buffer+writebuf.size, 0, superblock*record_size);
        }
        // Copy block into buffer
        memcpy(writebuf.buffer + writebuf.size, block, record_size);
        writebuf.size += (*maxbuf);
        
        // Normal particle block
        if (*ifl == 0)
            return;
        
        // Otherwise, this is the end of the shower. Dump the whole buffer,
        // zeroing out particle blocks for better compression if this is an
        // uninteresting shower.
        
        // Number of superblocks to write
        imax = ((writebuf.size/(*maxbuf))-1)/superblock + 1;
        for (i = 0; i < imax; i++) {
            for (j = 0; j < superblock; j++) {
                int offset = (i*superblock + j)*(*maxbuf);
                // If uninteresting, zero out everything except RUNH/RUNE, 
                // EVTH/EVTE 
                if (*zeroes &&
                    !(strncmp((char*)&writebuf.buffer[offset], "RUN", 3) == 0 ||
                    strncmp((char*)&writebuf.buffer[offset], "EVT", 3) == 0))
                    memset(&writebuf.buffer[offset], 0, record_size);
            }

            fwritempatap_(maxbuf, &superblock, &writebuf.buffer[i*superblock*(*maxbuf)]);
        }
        
        // Zero out everything copied over so far
        memset(writebuf.buffer, 0, imax*superblock*record_size);
        writebuf.size = 0;
    }
    
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

// ring buffer in C, where we have malloc()
// routines for ICECUBE1 FiFo by J. van Santen, U Wisconsin-Madison

typedef struct {
    size_t size, capacity;
    double *buffer, *begin, *end;
} ringbuffer_t ;

static ringbuffer_t corsika_ringbuffer = { 0, 0, NULL, NULL, NULL };

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
void ringbuffer_empty__(int *val)
#else
void ringbuffer_empty_(int *val)
#endif
{
    *val = (corsika_ringbuffer.size == 0);
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
void ringbuffer_clear__(int *maxlen)
#else
void ringbuffer_clear_(int *maxlen)
#endif
{
    ringbuffer_t *buf = &corsika_ringbuffer;
    // Shrink back down to a reasonable size
    if (buf->capacity > 2048*(*maxlen)) {
        buf->buffer = realloc(buf->buffer, 2048*(*maxlen)*sizeof(double));
        buf->capacity = 2048*(*maxlen);
    }
    buf->begin = buf->end = buf->buffer;
    buf->size = 0;
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#define max(a, b) (a > b ? a : b)

#if __REDHAT__
void ringbuffer_put__(double *particle, int *maxlen, int *at_end)
#else
void ringbuffer_put_(double *particle, int *maxlen, int *at_end)
#endif
{
    ringbuffer_t *buf = &corsika_ringbuffer;
    size_t record_size = (*maxlen)*sizeof(double);
    
    if (buf->size == buf->capacity) {
        size_t newsize = max(2048*(*maxlen), 2*buf->capacity);
        double *newbuf = malloc(newsize*sizeof(double));
        if (newbuf == NULL) {
            fprintf(stderr, "Could not allocate %zu bytes for ring buffer!\n", sizeof(double)*newsize);
         // exit(ENOMEM);
            exit(-1);
        }
        if (buf->buffer != NULL) {
            // NB: if buf->size == buf->capacity, then the buffer is full, i.e.
            // buf->end == buf->begin. The layout will look something like the
            // following sketch:
            //
            // |5 6 7 8 9 0 1 2 3 4|
            //  ^         ^         ^
            //  buf       end       buf+capacity
            //  -----------
            //    offset
            //
            // We want to copy the contents into a new, larger allocation in 
            // such a way that the empty region appears *after* the end pointer,
            // preserving the former ordering. Because addresses in the ring
            // buffer are periodic, this means that in general after the
            // expansion the empty region will appear between the *begin* and 
            // *end* pointers as in the following sketch:
            //
            //                     gap 
            //            _____________________
            // |5 6 7 8 9 * * * * * * * * * * * 0 1 2 3 4|
            //  ^         ^                     ^         ^
            //  buf       end                   begin     buf+capacity
            //  -----------
            //    offset

            // The distance between the start of the buffer and the end pointer 
            size_t offset = buf->end - buf->buffer;
            // The size of the new empty region
            size_t gap = newsize - buf->capacity;

            // Copy *offset* items to the front of the new buffer 
            memcpy(newbuf, buf->buffer, sizeof(double)*offset);
            buf->end = newbuf+offset;

            // Copy *capacity*-*offset* items to the back of the new buffer
            memcpy(newbuf+offset+gap, buf->buffer+offset,
                sizeof(double)*(buf->capacity-offset));
            buf->begin = newbuf+offset+gap;

            free(buf->buffer);
        } else {
            buf->begin = newbuf;
            buf->end   = newbuf;
        }
        
        buf->buffer = newbuf;
        buf->capacity = newsize;
    }
    
    if (*at_end == 1) {
        memcpy(buf->end, particle, record_size);
        buf->end += record_size;
        buf->size += record_size;

        // wrap around
        if (buf->end-buf->buffer == buf->capacity)
            buf->end = buf->buffer;
    } else {
        // NB: to write to the beginning of the ring buffer, just do
        // everything in reverse.
        if (buf->begin == buf->buffer)
            buf->begin = buf->buffer+buf->capacity;;

        buf->begin -= record_size;
        buf->size += record_size;
        memcpy(buf->begin, particle, record_size);
    }
}

#undef max

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#if __REDHAT__
void ringbuffer_get__(double *particle, int *maxlen)
#else
void ringbuffer_get_(double *particle, int *maxlen)
#endif
{
    ringbuffer_t *buf = &corsika_ringbuffer;
    size_t record_size = (*maxlen)*sizeof(double);
    
    if (buf->size == 0) {
        fprintf(stderr, "Nothing in the buffer!\n");
        exit(1);
        return;
    }
    
    memcpy(particle, buf->begin, record_size);
    buf->begin += record_size;
    buf->size -= record_size;
    
    // wrap around
    if (buf->begin-buf->buffer == buf->capacity)
        buf->begin = buf->buffer;
    
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#endif
