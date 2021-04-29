/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
NMUTILS.C - Functions for command line parsing, basic trace operations, 
            memory allocation, byte-swapping, and some others
*************************************************************************

Author: Nils Maercklin,
  RISSC, University of Naples, Italy, September 2008

Version: 2011-06-12

Modifications:
  2010-01-30 (NM): Moved many functions from "nmsaclib.c" to this file
  2010-03-16 (NM): Included functions operating on trace data
  2011-06-04 (NM): Included function to construct output file names

Notes:
  These functions are required by main NMSAC programs as well as by 
  some additional library functions, e.g. those in "nmsaclib.c".
  The function butterworth_bandpass() requires routines from "butterworth.c".
  All other functions use only standard C libraries.
*************************************************************************/

/* General include files */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#define STREQ(s,t)  (strcmp(s,t)     == 0)
#define STRCEQ(s,t) (strcasecmp(s,t) == 0)

/* Prototypes of functions for command line parsing */
int fgetpar(int argc, char **argv, char *par, float *val);
int igetpar(int argc, char **argv, char *par, int *val);
int sgetpar(int argc, char **argv, char *par, char **val);
int getflag(int argc, char **argv, char *par);

/* Prototypes of functions operating on trace data*/
void remove_mean(float *data, int nt);
void reverse_trace(float *data, int nt);
void butterworth_bandpass(int npoles, float flo, float fhi, float dt, \
    int nt, int zerophase, float *data);

/* Miscellaneous functions */
int out_file_name(char *fname, char *ext, char *dir, char **outfname);
void swab4(char *pt, int n);
float *fmalloc1(size_t n1);
float **fmalloc2(size_t n1, size_t n2);
void free2(void **p);
void error(char *fmt, ...);


/* Functions from "butterworth.c", required by butterworth_bandpass() */
void bfhighpass (int npoles, float f3db, int n, float p[], float q[]);
void bflowpass (int npoles, float f3db, int n, float p[], float q[]);


/************************************************************************/
/* Functions to parse command line arguments                            */
/************************************************************************/

int fgetpar(int argc, char **argv, char *par, float *val)
/************************************************************************
fgetpar - get a float value from the command line

Input:
argc    number of command line arguments
argv    array of command line arguments
par     command line argument name, e.g. "-a" (case-insensitive)

Output:
val     pointer to float variable (value of argv[i+1], if par==argv[i])
        returns position of argument name par in argv

Author: Nils Maercklin, 2008-09-05
*************************************************************************/
{
    register int i;

    for (i=1; i<argc-1; i++) {
        if (STRCEQ(argv[i], par)) {
            *val = atof(argv[i+1]);
            return i;
        }
    }
    return 0;
}



int igetpar(int argc, char **argv, char *par, int *val)
/************************************************************************
igetpar - get an integer value from the command line

Input:
argc    number of command line arguments
argv    array of command line arguments
par     command line argument name, e.g. "-a" (case-insensitive)

Output:
val     pointer to int variable (value of argv[i+1], if par==argv[i])
        returns position of argument name par in argv

Author: Nils Maercklin, 2008-09-05
*************************************************************************/
{
    register int i;

    for (i=1; i<argc-1; i++) {
        if (STRCEQ(argv[i], par)) {
            *val = atoi(argv[i+1]);
            return i;
        }
    }
    return 0;
}



int sgetpar(int argc, char **argv, char *par, char **val)
/************************************************************************
sgetpar - get a character string from the command line

Input:
argc    number of command line arguments
argv    array of command line arguments
par     command line argument name, e.g. "-a" (case-insensitive)

Output:
val     pointer to char array (argv[i+1], if par==argv[i])
        returns position of argument name par in argv

Author: Nils Maercklin, 2008-09-05
*************************************************************************/
{
    register int i;

    for (i=1; i<argc-1; i++) {
        if (STRCEQ(argv[i], par)) {
            val[0] = argv[i+1];
            return i;
        }
    }
    return 0;
}



int getflag(int argc, char **argv, char *par)
/************************************************************************
getflag - search for command line option and return its position

Input:
argc    number of command line arguments
argv    array of command line arguments
par     command line argument name, e.g. "-a" (case-insensitive)

Output:
val     returns position of argument name par in argv (0 if par not found)
*************************************************************************/
{
    register int i;

    for (i=1; i<argc; i++) {
        if (STRCEQ(argv[i], par)) {
            return i;
        }
    }
    return 0;
}



/************************************************************************/
/* Functions operating on trace data                                    */
/************************************************************************/

void remove_mean(float *data, int nt)
/************************************************************************
remove_mean - subtract the mean value from a data trace

Input:
data    array[nt] of data
nt      number of samples in data array

Output: modified data array

Author: Nils Maercklin, 2004-03-01
*************************************************************************/
{
    register int i;
    double mean=0.0;

    if (nt>0) {
        for (i=0; i<nt; i++) mean += (double)data[i];
        mean /= (double) nt;

        for (i=0; i<nt; i++) data[i] -= (float)mean;
    }
}



void reverse_trace(float *data, int nt)
/************************************************************************
reverse_trace - reverse trace in place

Input:
data    array[nt] of data
nt      number of samples in data array

Output: modified data array

Author: Nils Maercklin, 2004-03-01
*************************************************************************/
{
    register int i;
    register float tmp;

    for (i=0; i<nt/2; ++i) {
        tmp = data[i];
        data[i] = data[nt-1 - i];
        data[nt-1 - i] = tmp;
    }
}



void butterworth_bandpass(int npoles, float flo, float fhi, float dt, \
    int nt, int zerophase, float *data)
/************************************************************************
butterworth_bandpass - Butterworth bandpass filter

Input:
npoles  number of poles in Butterworth filters
flo     low-cut frequency in Hz  (0 < flo < F_Nyquist)
fhi     high-cut frequency in Hz (0 < fhi < F_Nyquist)
dt      time sampling interval in seconds
nt      number of samples in data array
zerophase  flag: 1 = zero-phase (two-pass) filter, 0 = minimum phase
data    array[nt] of input/output time-series

Output:
        modifies array data[]

Notes:
No filter is applied, if the normalized frequency is outside the valid range.
The actual filtering is done by routines from "butterworth.c", taken from 
the CWP/SU Seismic Unix package.

Author: Nils Maercklin, 2009-09-13
*************************************************************************/
{
    float f;   /* normalized frequency, e.g. 0.0 < flo*dt < 0.5) */

    /* Low-cut filter (high-pass) */
    f = flo * dt;
    if (npoles > 0 && f > 0.0 && f < 0.5) {
        bfhighpass(npoles, f, nt, data, data);
        if (zerophase) {
            reverse_trace(data, nt);
            bfhighpass(npoles, f, nt, data, data);
            reverse_trace(data, nt);
        }
    }

    /* High-cut filter (low-pass) */
    f = fhi * dt;
    if (npoles > 0 && f > 0.0 && f < 0.5) {
        bflowpass(npoles, f, nt, data, data);
        if (zerophase) {
            reverse_trace(data, nt);
            bflowpass(npoles, f, nt, data, data);
            reverse_trace(data, nt);
        }
    }
}



/************************************************************************/
/* Miscellaneous functions                                              */
/************************************************************************/

int out_file_name(char *fname, char *ext, char *dir, char **outfname)
/************************************************************************
out_file_name - set output file name based on input file and extension

Input:
fname    input file name
ext      file extension to be appended
dir      output directory (NULL = use directory from fname)

Output:
outfname output file name in the form "fname.ext" or "outdir/fname_base.ext"
         (pointer to character string char *)

Notes:
Remember to free outfname after usage. The function returns 0 on failure.

Author: Nils Maercklin, 2011-06-04
*************************************************************************/
{
    register int i;
    int len, ibase;

    /* File name with new output directory */
    ibase = 0;
    if (dir) {
        /* Get start of input name fname ("basename") */
        len = strlen(fname);
        for (i=len-1; i>=0; i--) {
            if (fname[i]=='/') {
                ibase = i + 1;
                break;
            }
        }

        /* Allocate space and begin constructing file name */
        len = strlen(dir) + strlen(fname+ibase) + strlen(ext) + 3;
        if (!(outfname[0]=malloc(len))) {
            return 0;
        }
        else {
            strcpy(outfname[0], dir);
            outfname[0]=strcat(outfname[0], "/");
            outfname[0]=strcat(outfname[0], fname+ibase);
        }
    }

    /* File name keeping input file directory */
    else {
        /* Allocate space and begin constructing file name */
        len = strlen(fname) + strlen(ext) + 3;
        if (!(outfname[0]=malloc(len))) {
            return 0;
        }
        else {
            strcpy(outfname[0], fname);
        }
    }

    /* Append file extension */
    outfname[0]=strcat(outfname[0], ".");
    outfname[0]=strcat(outfname[0], ext);

    return 1;
}




void swab4(char *pt, int n) 
/************************************************************************
swab4 - reverse byte order for 4-byte float/integer

Input:
pt      pointer to byte array
n       number of bytes

Output: n swapped bytes at *pt

Author: Lupei Zhu, 1996-12-03 (function taken from "sacio.c")
*************************************************************************/
{
    register int i;
    char temp;

    for(i=0; i<n; i+=4) {
        temp    = pt[i+3];
        pt[i+3] = pt[i];
        pt[i]   = temp;
        temp    = pt[i+2];
        pt[i+2] = pt[i+1];
        pt[i+1] = temp;
    }
}




float *fmalloc1(size_t n1)
/************************************************************************
fmalloc1 - allocate 1-D array of floats and initialize with zeros

Input:
n1      size of arrays (number of floats)

Output: returns pointer to allocated array (NULL on failure)
*************************************************************************/
{
    float *p=NULL;

    if ((p=malloc(n1*sizeof(float)))==NULL) {
        return NULL;
    }
    else {
        memset((void *)p, 0, n1*sizeof(float));
        return p;
    }
}




void **alloc2(size_t n1, size_t n2, size_t size);
float **fmalloc2(size_t n1, size_t n2)
/************************************************************************
fmalloc2 - allocate 2-D array of floats and initialize with zeros

Input:
n1, n2   size of array (number of floats)

Output: returns pointer to allocated array (NULL on failure)
*************************************************************************/
{
    float **p=NULL;

    if ((p=(float**)alloc2(n1, n2, sizeof(float)))==NULL) {
        return NULL;
    }
    else {
        memset((void *)p[0], 0, n1*n2*sizeof(float));
        return p;
    }
}

void **alloc2(size_t n1, size_t n2, size_t size)
{
    size_t i2;
    void **p=NULL;

    if ((p=(void**)malloc(n2*sizeof(void*)))==NULL) 
        return NULL;

    if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
        free(p);
        return NULL;
    }

    for (i2=0; i2<n2; i2++)
        p[i2] = (char*)p[0]+size*n1*i2;

    return p;
}

void free2(void **p)
{
    free(p[0]);
    free(p);
}




void error(char *fmt, ...)
/************************************************************************
error - write an error message and exit program (wrapper for fprintf())
*************************************************************************/
{
    va_list args;

    fprintf(stderr, "\n");
    va_start(args,fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    exit(EXIT_FAILURE);
}

/* END OF FILE */
