/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
NMSACLIB.H - Definitions of SAC header and file utilities ("nmsaclib.c"), 
             prototypes of library functions used by NMSAC programs
*************************************************************************

Author: Nils Maercklin,
  RISSC, University of Naples, Italy, September 2008

Version: 2011-06-12

Modifications:
  2010-01-30 (NM): Changed SAC header definition file, and included
                   prototypes of other library functions
  2010-03-16 (NM): Updated prototypes of other library functions

Notes:
  For a description of defined functions, their usage, dependencies, 
  and of possible limitations see "nmsaclib.c" or other files listed
  below. The SAC header is defined in "nmsac.h", but all SAC-related
  functions defined below are also compatible with the SAC header 
  definitions in "utils/sac.h" distributed with SAC 101.2 (2008).
  Routines or definitions from SAC itself are not required here.

  For source files of main NMSAC programs it should be sufficient
  to include this header file "nmsaclib.h".
  

Reference (SAC file format):
  SAC - Seismic Analysis Code, User's manual, 
  http://www.iris.edu/manuals/sac/manual.html (last checked: 2008-02-27)
*************************************************************************/



/************************************************************************/
/* Include files and definitions                                        */
/************************************************************************/

/* General include files */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <float.h>
#include <math.h>

/* Include SAC header definition file "nmsac.h" */
#include "nmsac.h"

/* Constants */
#define SAC_HEADER_NVHDR 6


/* Check type of SAC header fields */
#define IS_SACHDR_NUMERIC(id) \
    (( (id)>=0 && (id)<=SAC_HEADER_LOGICAL_MAX ) ? 1 : 0 )
#define IS_SACHDR_FLOAT(id) \
    (( (id)>=SAC_HEADER_FLOAT_MIN && (id)<=SAC_HEADER_FLOAT_MAX ) ? 1 : 0 )
#define IS_SACHDR_INT(id) \
    (( (id)>=SAC_HEADER_INT_MIN && (id)<=SAC_HEADER_LOGICAL_MAX ) ? 1 : 0 )
#define IS_SACHDR_CHAR(id) \
    (( (id)>=SAC_HEADER_CHAR_MIN && (id)<=SAC_HEADER_CHAR_MAX ) ? 1 : 0 )

#define IS_SACHDR_TMARK(id) (( (id)==8 || ( (id)>=SAC_HEADER_TMARK_POSITION \
    && (id)<=SAC_HEADER_TMARK_POSITION+10 )) ? 1 : 0 )


/* Check for byte order of SAC file (based on header version number) */
#define NEED_BYTESWAP(i) ( ( ((i)==SAC_HEADER_INT_UNDEFINED) \
    || ((i)>=0 && (i)<= SAC_HEADER_NVHDR) ) ? 0 : 1 )


/* Abbreviations */
#ifndef PI
#define PI          ( 3.141592653589793 )
#endif
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#ifndef SGN
#define SGN(x) (((x) < 0.0) ? -1.0 : 1.0)
#endif
#define NINT(x)     ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define STREQ(s,t)  (strcmp(s,t)     == 0)
#define STRCEQ(s,t) (strcasecmp(s,t) == 0)



/************************************************************************/
/* Prototypes of functions in "nmsaclib.c"                              */
/************************************************************************/

/* Prototypes of functions related to SAC headers */
int sac_hdr_index(const char *key);
int sac_hdr_offset(int index);
int sac_hdr_enum(const char *name);
char *sac_enum_name(int value, int desc);
char *sac_hdr_name(int index);

float get_sac_hdr_float(const SACHEAD *hd, int index);
int get_sac_hdr_int(const SACHEAD *hd, int index);
int get_sac_hdr_string(const SACHEAD *hd, int index, int rb, char *string);

void put_sac_hdr_float(SACHEAD *hd, int index, const float val);
void put_sac_hdr_int(SACHEAD *hd, int index, const int val);
void put_sac_hdr_string(SACHEAD *hd, int index, const char *string);

void set_sac_depminmax(SACHEAD *hd, const float *data);


/* Prototypes of functions for SAC file/trace input/output */
int read_sacbin(FILE *fp, SACHEAD *hd, float **data, int *swap);
int read_sacbin_file(char *filename, SACHEAD *hd, float **data, int *swap);
int write_sacbin(FILE *fp, SACHEAD hd, float *data, int swap);
int write_sacbin_file(char *filename, SACHEAD hd, float *data, int swap);



/************************************************************************/
/* Prototypes of functions in "nmutils.c"                               */
/************************************************************************/

/* Prototypes of functions for command line parsing */
int fgetpar(int argc, char **argv, char *par, float *val);
int igetpar(int argc, char **argv, char *par, int *val);
int sgetpar(int argc, char **argv, char *par, char **val);
int getflag(int argc, char **argv, char *par);

/* Prototypes of functions for basic trace operations */
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



/************************************************************************/
/* Prototypes of functions in "nmpolar.c"                               */
/************************************************************************/

/* Prototypes of functions related to polarization analyses */
void eig_jacobi(float **a, float d[], float **v, int n);
void sort_eigenvalues(float d[], float **v, int n);

float covar(float *data1, float *data2, int istart, int iwl);
void do_eigen(float **data, int nt, int iwl, float **ddata, float **vdata);
void do_eigen_zm(float **data, int nt, int iwl, float **ddata, float **vdata);

float calc_rl(float *d, float rlq, int opt);
float calc_ellip(float *d, int i1, int i2);
float calc_f1(float *d);
float calc_l1(float *d);
float calc_plan(float *d);
float calc_tau(float *d);
float calc_er(float *d);
float calc_norminc(float *v);
float calc_theta(float *v, int opt);
float calc_phi(float *v, int opt);



/************************************************************************/
/* Prototypes of functions in "nmetime.c"                               */
/************************************************************************/

/* Prototypes of functions for epoch time conversions */
int ndays_yearmonth(int year, int month);
int hdate_to_jday(int year, int month, int day);
int jday_to_hdate(int year, int jday, int *month, int *day);
double date_to_etime(int eyear, int year, int jday, int hour, int minute, \
    float sec);
void etime_to_date(double etime, int eyear, int *year, int *jday, \
    int *hour, int *minute, float *sec);



/************************************************************************/
/* Prototypes of functions in "nmgeo.c"                                 */
/************************************************************************/

/* Function prototypes for geographical calculations */
double gc_azimuth(double lat1, double lon1, double lat2, double lon2);
double gc_delta(double lat1, double lon1, double lat2, double lon2);
void gc_inter(double *olat, double *olon, \
    double lat1, double lon1, double lat2, double lon2, double frac);
double gc_dist(double lat1, double lon1, double lat2, double lon2,
    double a, double f);



/************************************************************************/
/* Prototypes of functions in "butterworth.c" (from CWP/SU 41 (2008)    */
/************************************************************************/

/* Prototypes for design and application of Butterworth filters */
void bfhighpass (int npoles, float f3db, int n, float p[], float q[]);
void bflowpass (int npoles, float f3db, int n, float p[], float q[]);
void bfdesign (float fpass, float apass, float fstop, float astop,
    int *npoles, float *f3db);

/* END OF FILE */
