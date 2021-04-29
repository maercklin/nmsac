/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
NMSACLIB.C - SAC header and file utilities
*************************************************************************

Author: Nils Maercklin,
  RISSC, University of Naples, Italy, September 2008

Version: 2011-06-12

Modifications:
  2008-10-05 (NM): Code cleanup and minor other modifications
  2009-09-09 (NM): Added set_sac_depminmax(), and fmalloc2() etc.
  2010-01-30 (NM): Added sac_hdr_enum() etc., changed SAC header file, 
                   and moved not-SAC-specific functions to "nmutils.c"

Notes:
  Only evenly-sampled, binary SAC waveform files (LEVEN=1) are fully 
  supported, and the SAC header version should be NVHDR=6.

  The SAC header definition is given in "nmsac.h". Some global variables
  are defined in "nmsacnam.h", which must not be included in other source
  files. All functions defined below are also compatible with the 
  SAC header file "utils/sac.h" distributed with SAC 101.2 (2008).
  Routines from SAC itself are not required by these functions.

Reference (SAC file format):
  SAC - Seismic Analysis Code, User's manual, 
  http://www.iris.edu/manuals/sac/manual.html (last checked: 2008-02-27)
*************************************************************************/



/************************************************************************/
/* Include files, definitions, and function prototypes                  */
/************************************************************************/

/* General include files */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <float.h>
#include <math.h>

/* SAC-related definitions */
#include "nmsac.h"     /* SAC header definition */
#include "nmsacnam.h"  /* Global variables, header field names and values */
#include "nmutils.h"   /* Utility functions */
#include "nmsaclib.h"  /* SAC-related functions */


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
/* Functions related to SAC binary trace headers                        */
/************************************************************************/

int sac_hdr_index(const char *key)
/************************************************************************
sac_hdr_index - get the index of a SAC header field

Input:
key     SAC header field name, see SACHEAD typedef for details

Output: returns field index (or -1 on failure)

Author: Nils Maercklin, 2006-10-26 (modified for "sac.h", 2008-09-03)
*************************************************************************/
{
    register int i;

    for (i = 0; i <= SAC_HEADER_FIELDS; i++) {
        if (STRCEQ(SacHeaderName[i], key)) return i;
    }

    /* key not found (return -1) */
    return -1;
}



int sac_hdr_offset(int index)
/************************************************************************
sac_hdr_offset - get the position (byte offset) of a SAC header field

Input:
index   SAC header field index (between 0 and SAC_HEADER_FIELDS)

Output: returns byte offset

Author: Nils Maercklin, 2008-09-03
*************************************************************************/
{
    int offs = 0;

    if (index < 0 || index > SAC_HEADER_FIELDS) {
        error("%s: invalid header field index %d\n", __func__, index);
    }

    if (index < SAC_HEADER_CHAR_MIN) {
        offs = 4*index;
    }
    else {
        offs = SAC_HEADER_SIZE_NUMBERS + 8*(index - SAC_HEADER_CHAR_MIN);
    }

    return offs;
}



int sac_hdr_enum(const char *name)
/************************************************************************
sac_hdr_enum - get the numerical value of a SAC header enum name

Input:
name    SAC header enum name, see SacHeaderEnums in "nmsacnam.h" or "sac.h"

Output: returns enum value (or -1 on failure)

Author: Nils Maercklin, 2008-10-05
*************************************************************************/
{
    register int i;

    for (i = 0; i <= SacHeaderEnumsLength; i++) {
        if (STRCEQ(SacHeaderEnums[i], name)) return i;
    }

    /* enum name not found (return -1) */
    return -1;
}



char *sac_enum_name(int value, int desc)
/************************************************************************
sac_enum_name - get the name of a SAC header enum field from its value

Input:
value   SAC header enum value, see "nmsacnam.h" or "sac.h" for details
desc    enum name flag, 0 = short name, 1 = long description

Output: returns enum name (or dummy string on failure)

Author: Nils Maercklin, 2010-01-30
*************************************************************************/
{
    if (desc) {
        if (value < 0 || value > SacHeaderEnumsLength) {
            return "Unknown ENUM";
        }
        else if (value == SAC_HEADER_INT_UNDEFINED) {
            return "Undefined";            
        }
        else {
            return SacHeaderEnumsDescription[value];
        }
    }
    else {
        if (value < 0 || value > SacHeaderEnumsLength) {
            return "ENUM?";
        }
        else if (value == SAC_HEADER_INT_UNDEFINED) {
            return "UNDEF";            
        }
        else {
            return SacHeaderEnums[value];
        }
    }
}



char *sac_hdr_name(int index)
/************************************************************************
sac_hdr_name - get the name of a SAC header field from its index

Input:
index   SAC header field index, see SACHEAD typedef for details

Output: returns field name (or "null" on failure)

Author: Nils Maercklin, 2010-01-30
*************************************************************************/
{
    if (index < 0 || index > SAC_HEADER_FIELDS) {
        return "null";
    }
    else {
        return SacHeaderName[index];
    }
}



float get_sac_hdr_float(const SACHEAD *hd, int index)
/************************************************************************
get_sac_hdr_float - get float value from SAC binary header

Input:
hd      SAC header structure   (see SACHEAD typedef)
index   SAC header field index (between 0 and SAC_HEADER_FIELDS)

Output: returns header value cast to float

Author: Nils Maercklin, 2008-09-07
*************************************************************************/
{
    int offs = sac_hdr_offset(index);
    char *tp = (char *) hd;
    char buff[9];

    if (index < 0 || index > SAC_HEADER_FIELDS) {
        error("%s: invalid header field index %d\n", __func__, index);
    }

    if (IS_SACHDR_FLOAT(index)) {
        return *((float*) (tp+offs));
    }
    else if (IS_SACHDR_INT(index)) {
        return (float) *((int*) (tp+offs));
    }
    else {
        strncpy(buff, tp+offs, 8);
        buff[8]='\0';
        return atof(buff);
    }
}



int get_sac_hdr_int(const SACHEAD *hd, int index)
/************************************************************************
get_sac_hdr_int - get integer value from SAC binary header

Input:
hd      SAC header structure   (see SACHEAD typedef)
index   SAC header field index (between 0 and SAC_HEADER_FIELDS)

Output: returns header value cast to int

Author: Nils Maercklin, 2008-09-07
*************************************************************************/
{
    int offs = sac_hdr_offset(index);
    char *tp = (char *) hd;
    char buff[9];

    if (index < 0 || index > SAC_HEADER_FIELDS) {
        error("%s: invalid header field index %d\n", __func__, index);
    }

    if (IS_SACHDR_FLOAT(index)) {
        return (int) *((float*) (tp+offs));
    }
    else if (IS_SACHDR_INT(index)) {
        return *((int*) (tp+offs));
    }
    else {
        strncpy(buff, tp+offs, 8);
        buff[8]='\0';
        return atoi(buff);
    }
}



int get_sac_hdr_string(const SACHEAD *hd, int index, int rb, char *string)
/************************************************************************
get_sac_hdr_string - get character string from SAC binary header

Input:
hd      SAC header structure   (see SACHEAD typedef)
index   SAC header field index (between 0 and SAC_HEADER_FIELDS)
rb      remove-blanks flag: 0 = keep string unmodified, 1 = remove
        trailing blanks, 2 = remove trailing and leading blanks, 
string  character string buffer (size larger than 16 characters)

Output: returns string length
string  character string from header, including "\0"

Author: Nils Maercklin, 2008-10-05
*************************************************************************/
{
    register int i;
    int offs = sac_hdr_offset(index);
    int slen = 0;
    char *tp = (char *) hd;
    char buff[30];

    if (index < 0 || index > SAC_HEADER_FIELDS) {
        error("%s: invalid header field index %d\n", __func__, index);
    }

    if (IS_SACHDR_FLOAT(index)) {
        slen = 16;
        (void)sprintf(buff, "%-16f", get_sac_hdr_float(hd,index));
    }
    else if (IS_SACHDR_INT(index)) {
        slen = 16;
        (void)sprintf(buff, "%-16d", get_sac_hdr_int(hd,index));
    }
    else {
        slen = (index==SAC_HEADER_CHAR_DOUBLE) ? 16 : 8;
        strncpy(buff, tp+offs, slen);
        buff[slen+1]='\0';
    }

    strncpy(string, buff, (slen<16)?slen:16);
    string[slen]='\0';

    /* Remove trailing and/or leading blanks */
    if (rb) {
        for (i=slen-2;i>0;i--) {
            if (string[i]==' ') string[i]='\0';
            else break;
        }
    }
    if (rb>1) {
        while (string[0]==' ') {
            for (i=0;i<slen;i++) string[i]=string[i+1];
        }
    }

    return strlen(string);
}



void put_sac_hdr_float(SACHEAD *hd, int index, const float val)
/************************************************************************
put_sac_hdr_float - write a float value to SAC binary header

Input:
hd      SAC header structure   (see SACHEAD typedef)
index   SAC header field index (between 0 and SAC_HEADER_FIELDS)
val     float value to be written to header field specified by index

Output:
hd      modified SAC header

Author: Nils Maercklin, 2008-09-07
*************************************************************************/
{
    int offs = sac_hdr_offset(index);
    char *tp = (char *) hd;
    char buff[30];

    if (IS_SACHDR_FLOAT(index)) {
        *((float*) (tp+offs)) = val;
    }
    else if (IS_SACHDR_INT(index)) {
        *((int*) (tp+offs)) = NINT(val);
    }
    else {
        (void)sprintf(buff, "%-8f", val);
        strncpy(tp+offs, buff, 8);
    }
}



void put_sac_hdr_int(SACHEAD *hd, int index, const int val)
/************************************************************************
put_sac_hdr_int - write an integer value to SAC binary header

Input:
hd      SAC header structure   (see SACHEAD typedef)
index   SAC header field index (between 0 and SAC_HEADER_FIELDS)
val     int value to be written to header field specified by index

Output:
hd      modified SAC header

Author: Nils Maercklin, 2008-09-07
*************************************************************************/
{
    int offs = sac_hdr_offset(index);
    char *tp = (char *) hd;
    char buff[30];

    if (IS_SACHDR_FLOAT(index)) {
        *((float*) (tp+offs)) = (float) val;
    }
    else if (IS_SACHDR_INT(index)) {
        *((int*) (tp+offs)) = val;
    }
    else {
        (void)sprintf(buff, "%-8d", val);
        strncpy(tp+offs, buff, 8);
    }
}



void put_sac_hdr_string(SACHEAD *hd, int index, const char *string)
/************************************************************************
put_sac_hdr_string - write a character string to SAC binary header

Input:
hd      SAC header structure   (see SACHEAD typedef)
index   SAC header field index (between 0 and SAC_HEADER_FIELDS)
string  character string to be written to header field specified by index

Output:
hd      modified SAC header

Author: Nils Maercklin, 2008-10-05
*************************************************************************/
{
    int offs = sac_hdr_offset(index);
    char *tp = (char *) hd;
    int slen = 0;
    int ilen = 0;
    int i    = 0;

    if (IS_SACHDR_FLOAT(index)) {
        *((float*) (tp+offs)) = atof(string);
    }
    else if (IS_SACHDR_INT(index)) {
        *((int*) (tp+offs)) = atoi(string);
    }
    else {
        slen = (index==SAC_HEADER_CHAR_DOUBLE) ? 16 : 8;
        ilen = strlen(string);
        if (ilen>=slen) {
            strncpy(tp+offs, string, slen);
        }
        else {
            strncpy(tp+offs, string, slen);
            for (i=ilen; i<slen; i++) *(tp+offs+i) = ' ';
        }
    }
}



void set_sac_depminmax(SACHEAD *hd, const float *data)
/************************************************************************
set_sac_depminmax - compute and set SAC header demin/depmax/depmean

Input:
hd      SAC header structure   (see SACHEAD typedef)
data    data trace, pointer to array of floats

Output:
hd      modified SAC header

Author: Nils Maercklin, 2008-10-31
*************************************************************************/
{
    float  fmin = 0.0;
    float  fmax = 0.0;
    double dmean= 0.0;
    int    npts = 0;
    register int i;

    npts = hd->npts;

    if (npts==SAC_HEADER_INT_UNDEFINED) {
        hd->depmin = SAC_HEADER_FLOAT_UNDEFINED;
        hd->depmax = SAC_HEADER_FLOAT_UNDEFINED;
        hd->depmen = SAC_HEADER_FLOAT_UNDEFINED;
    }
    else {
        for (i=0; i<npts; i++) {
            dmean += (double) data[i];

            if (!i || data[i] > fmax) fmax = data[i];
            if (!i || data[i] < fmin) fmin = data[i];
        }
        if (npts) dmean /= (double) npts;

        hd->depmin = fmin;
        hd->depmax = fmax;
        hd->depmen = (float) dmean;
    }
}



/************************************************************************/
/* Functions to read and write evenly-sampled SAC binary traces         */
/************************************************************************/

int read_sacbin(FILE *fp, SACHEAD *hd, float **data, int *swap)
/************************************************************************
read_sacbin - read SAC binary trace from input stream

Input:
fp      FILE pointer, input stream
hd      SAC header structure (see SACHEAD typedef)
data    pointer to array of floats

Output:
swap    byte-swapping flag, 1 if byte order reversed, else 0
        returns 1 if successful, 0 on failure

Note:   Memory allocation for *data within this function (hd.npts floats);
        remember to free the space, if *data is no longer needed.

Author: Nils Maercklin, 2008-09-07
*************************************************************************/
{
    *swap=0;

    /* Read SAC header */
    if (fread(hd, sizeof(SACHEAD), 1, fp) != 1) return 0;

    /* Swap bytes, if necessary */
    if (NEED_BYTESWAP(hd->nvhdr)) {
        swab4((char*) hd, SAC_HEADER_SIZE_NUMBERS);
        *swap=1;
    }
   
    /* Allocate space for SAC trace data */
    if (!(data[0]=fmalloc1(hd->npts))) {
        fprintf(stderr, "%s: can't allocate space, npts=%d\n", \
            __func__, hd->npts);
        return 0;
    }

    /* Read SAC trace data */
    if (fread(data[0], sizeof(float), hd->npts, fp) != hd->npts) {
        free(data[0]);
        return 0;
    }

    /* Swap bytes */
    if ((*swap)) swab4((char*) data[0], (hd->npts)*sizeof(float));

    return 1;
}



int read_sacbin_file(char *filename, SACHEAD *hd, float **data, int *swap)
/************************************************************************
read_sacbin_file - read an evenly-sampled SAC binary waveform file

Input:
filename  name of a SAC binary waveform file
hd      SAC header structure (see SACHEAD typedef)
data    pointer to array of floats

Output:
swap    byte-swapping flag, 1 if byte order reversed, else 0
        returns 1 if successful, 0 on failure

Note:   Memory allocation for *data within this function (hd.npts floats);
        remember to free the space, if *data is no longer needed.

Author: Nils Maercklin, 2008-09-07
*************************************************************************/
{
    FILE *fp=NULL;

    /* Open input file */
    if (!(fp=fopen(filename, "r"))) return 0;

    /* Write SAC binary trace */
    if (!read_sacbin(fp, hd, data, swap)) {
        fclose(fp);
        return 0;
    }

    /* Close output file */
    fclose(fp);

    return 1;
}



int write_sacbin(FILE *fp, SACHEAD hd, float *data, int swap)
/************************************************************************
write_sacbin - write SAC binary trace to output stream

Input:
fp      FILE pointer, output stream
hd      SAC header structure (see SACHEAD typedef)
data    array of hd.npts floats
swap    byte-swapping flag, 1 = reverse bytes, 0 = do not

Output: returns 1 if successful, 0 on failure

Author: Nils Maercklin, 2008-09-07
*************************************************************************/
{
    int npts=hd.npts;

    /* Byte-swapping */
    if (swap) {
        swab4((char*) &hd, SAC_HEADER_SIZE_NUMBERS);
        swab4((char*) data, npts*sizeof(float));
    }

    /* Write SAC header */
    if (fwrite(&hd, sizeof(SACHEAD), 1, fp) != 1) return 0;

    /* Write SAC trace data */
    if (fwrite(data, sizeof(float), npts, fp) != npts) return 0;

    return 1;
}



int write_sacbin_file(char *filename, SACHEAD hd, float *data, int swap)
/************************************************************************
write_sacbin_file - write an evenly-sampled SAC binary waveform file

Input:
filename  name of output SAC binary waveform file
hd      SAC header structure (see SACHEAD typedef)
data    array of hd.npts floats
swap    byte-swapping flag, 1 = reverse bytes, 0 = do not

Output: returns 1 if successful, 0 on failure

Author: Nils Maercklin, 2008-09-07
*************************************************************************/
{
    FILE *fp=NULL;

    /* Open output file */
    if (!(fp=fopen(filename, "w"))) return 0;

    /* Write SAC binary trace */
    if (!write_sacbin(fp, hd, data, swap)) {
        fclose(fp);
        return 0;
    }

    /* Close output file */
    fclose(fp);

    return 1;
}

/* END OF FILE */
