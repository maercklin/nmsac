/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
NMUTILS.H - Definitions of utility functions in "nmutils.c"
*************************************************************************

Author: Nils Maercklin,
  RISSC, University of Naples, Italy, September 2008

Version: 2011-06-12

Modifications:
  2010-01-30 (NM): Moved many definitions from "nmsaclib.c" to this file
  2010-03-16 (NM): Included functions operating on trace data
  2011-06-04 (NM): Included function to construct output file names

Notes:
  For details on defined functions refer to the file "nmutils.c".

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

/* END OF FILE */
