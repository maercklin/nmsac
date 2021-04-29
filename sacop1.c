/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACOP1 - Unary arithmetic Operations on SAC trace data
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, September 2009

Version: 2011-06-12

Modifications:
    2009-09-15 (NM): First version of this code
    2011-06-12 (NM): Added user-specified output directory

Notes: Input files must be evenly-sampled time series files in SAC 
    binary format (assuming NVHDR=6, IFTYPE=ITIME (1), LEVEN=1).
*************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "nmsaclib.h"     /* SAC-related functions, and other stuff */


/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACOP1 - Unary arithmetic Operations on SAC trace data                ",
"                                                                       ",
" Usage: sacop1 [-p] operation  [parameters] -f sac_files               ",
"        sacop1 <stdin [-p] operation [parameters] > stdout             ",
"                                                                       ",
" Parameters and defaults:                                              ",
"   -f        list of SAC binary waveform files (or read from stdin)    ",
"   -p (none) one of the arithmetic operations listed below             ",
"   -s   1.0  trace scale factor applied before operation given by -p   ",
"   -d        flag: input in degrees instead of radians (sin/cos/tan)   ",
"   -o        flag: overwrite input files (default appends suffix)      ",
"   -dir      output directory (default is same as input directory)     ",
"                                                                       ",
" Arithmetic operations:                                                ",
"   abs, amp  absolute amplitude of trace (modulus)                     ",
"   sqr, ssqr squared amplitude and signed squared amplitude            ",
"   ssqrt     signed square root                                        ",
"   slog10    signed common logarithm                                   ",
"   sin, cos  trigonometric functions sine and cosine (input in rad)    ",
"   tan       trigonometric function tangens (input in rad)             ",
"   inv       inverse of trace amplitudes                               ",
"   rmean     remove mean value                                         ",
"   norm      normalize amplitudes by absolute maximum                  ",
"   sum       running sum trace integration                             ",
"   diff      running difference trace differentiation                  ",
"   nop       no operation                                              ",
"                                                                       ",
" One operation must be specified. Use e.g. nop for trace scaling only. ",
" By default, the specified operation name is appended to all output    ",
" files. If the overwrite flag -o is set, original data are lost.       ",
"                                                                       ",
" NM, 2011-06-12                                                        ",
"                                                                       ",
NULL};

/* Prototypes of functions used internally */
void normalize_trace(float *data, int nt);
void differentiate(float *data, int nt, float dt);

#define RAD2DEG  57.29577951
#define DEG2RAD  0.017453292

int 
main(int argc, char **argv)
{
    /* Variables */
    register int iarg, it; /* loop indices */
    SACHEAD hd;          /* SAC header */
    float *data=NULL;    /* SAC input trace data */
    char *outfname=NULL; /* output file name */
    char *outdir=NULL;   /* optional output directory name */
    char *op=NULL;       /* arithmetic operation name */
    int isfile=0;        /* file names flag */
    int swap;            /* byte-swapping flag */
    int overwrite=0;     /* overwrite flag */
    int degrees=0;       /* zero-phase flag */
    float a;             /* internal scale factor or temporary float */
    float s;             /* user-specified scale factor */


    /* Print documentation */
    if (argc==1) {
        int i=0;
        while (sdoc[i]) fprintf(stderr, "%s\n", sdoc[i++]);
        return (EXIT_FAILURE);
    }

    /* Command line parameters */
    if (!sgetpar(argc, argv, "-p", &op))         op=argv[1];
    if (!fgetpar(argc, argv, "-s", &s))          s=1.0;
    if (!(overwrite=getflag(argc, argv, "-o")))  overwrite=0;
    if (!(degrees=getflag(argc, argv, "-d")))    degrees=0;
    if (!(isfile=getflag(argc, argv, "-f")))     isfile=0;
    if (!sgetpar(argc, argv, "-dir", &outdir))   outdir=NULL;

    /* Check parameters */
    if (!isfile && isatty(1)) \
        error("%s: can't write a SAC binary file to tty\n", argv[0]);

    /* Loop over SAC files */
    for (iarg=isfile+1; (iarg<argc && argv[iarg][0]!='-') || !isfile; iarg++) {

        /* Read SAC binary file */
        if (isfile) {
            if (!read_sacbin_file(argv[iarg], &hd, &data, &swap)) \
                error("%s: can't read %s\n", argv[0], argv[iarg]);
        }
        else {
            if (!read_sacbin(stdin, &hd, &data, &swap)) break;
        }

        /* Scale trace */
        if (s!=1.0) {
            for (it=0; it<hd.npts; it++) data[it] *= s;
        }

        /* Apply arithmetic operation */
        if (STRCEQ(op, "abs") || STRCEQ(op, "amp")) {
            for (it=0; it<hd.npts; it++) {
                data[it] = fabs(data[it]);
            }
        }
        else if (STRCEQ(op, "sqr")) {
            for (it=0; it<hd.npts; it++) {
                data[it] = data[it]*data[it];
            }
        }
        else if (STRCEQ(op, "ssqr")) {
            for (it=0; it<hd.npts; it++) {
                data[it] = SGN(data[it]) * (data[it]*data[it]);
            }
        }
        else if (STRCEQ(op, "ssqrt")) {
            for (it=0; it<hd.npts; it++) {
                data[it] = SGN(data[it]) * sqrt( fabs(data[it]) );
            }
        }
        else if (STRCEQ(op, "sin")) {
            a = (degrees) ? DEG2RAD : 1.0;
            for (it=0; it<hd.npts; it++) {
                data[it] = sin( a*data[it] );
            }
        }
        else if (STRCEQ(op, "cos")) {
            a = (degrees) ? DEG2RAD : 1.0;
            for (it=0; it<hd.npts; it++) {
                data[it] = cos( a*data[it] );
            }
        }
        else if (STRCEQ(op, "tan")) {
            a = (degrees) ? DEG2RAD : 1.0;
            for (it=0; it<hd.npts; it++) {
                data[it] = tan( a*data[it] );
            }
        }
        else if (STRCEQ(op, "inv")) {
            for (it=0; it<hd.npts; it++) {
                data[it] = (data[it]) ? 1.0/data[it] : 0.0;
            }
        }
        else if (STRCEQ(op, "slog10")) {
            for (it=0; it<hd.npts; it++) {
                a = fabs(data[it]);
                if (a>0.0) data[it] = SGN(data[it]) * log10(a);
                else       data[it] = 0.0;
            }
        }
        else if (STRCEQ(op, "norm")) {
            normalize_trace(data, hd.npts);
        }
        else if (STRCEQ(op, "rmean") || STRCEQ(op, "avg")) {
            remove_mean(data, hd.npts);
        }
        else if (STRCEQ(op, "diff")) {
            a = hd.delta;
            if (a==SAC_HEADER_FLOAT_UNDEFINED || a==0.0) a = 1.0;
            differentiate(data, hd.npts, a);
        }
        else if (STRCEQ(op, "sum") || STRCEQ(op, "int")) {
            a = hd.delta;
            if (a==SAC_HEADER_FLOAT_UNDEFINED || a==0.0) a = 1.0;
            data[0] *= a;
            for (it=1; it<hd.npts; it++) {
                data[it] *= a;
                data[it] += data[it-1];
            }
        }
        else if (STRCEQ(op, "nop")) {
            ;
        }
        else {
            free(data);
            error("%s: unknown operation \"%s\", doing nothing\n",\
                argv[0], op);
        }

        /* Re-compute depmen/depmin/depmax */
        set_sac_depminmax(&hd, data);

        /* Write SAC binary files... */
        if (isfile) {
            if (overwrite) {
                outfname = argv[iarg];
            }
            else {
                if (!out_file_name(argv[iarg], op, outdir, &outfname)) \
                    error("%s: can't set output file name\n", argv[0]);
            }
            if (!write_sacbin_file(outfname, hd, data, swap))\
                error("%s: can't write %s\n", argv[0], outfname);
            if (!overwrite) free(outfname);
        }
        /* ... or write SAC binary to stdout */
        else {
            if (!write_sacbin(stdout, hd, data, swap)) \
                error("%s: can't write SAC trace to stdout\n", argv[0]);
        }

        /* Free space */
        free(data);

    }

    return EXIT_SUCCESS;
}

/************************************************************************/
/* Function used internally                                             */
/************************************************************************/

void normalize_trace(float *data, int nt)
/************************************************************************
normalize_trace - normalize trace data by maximum value

Input:
data    array[nt] of data
nt      number of samples in data array

Output: modified data array

Author: Nils Maercklin, 2000 - 2009
*************************************************************************/
{
    register int i;
    float x, max;

    max = 0.0;
    for (i=0; i<nt; i++) {
        x = fabs(data[i]);
        if (max < x) max = x;
    }

    if (max>0.0) {
        for (i=0; i<nt; i++) data[i] /= max;
    }
}



void differentiate(float *data, int nt, float dt)
/************************************************************************
differentiate - running difference trace differentiation

Input:
data    array[nt] of data
nt      number of samples in data array
dt      sampling interval

Output: modified data array

Author: Nils Maercklin, 2009-09-15  (Credit: CWP/SU)
*************************************************************************/
{
    register int i;
    float *tmpdata=NULL;

    /* Allocate space and copy data */
    tmpdata = fmalloc1(nt);
    memcpy((void *)tmpdata, (const void *) data, nt*sizeof(float));

    if (dt==0.0) dt = 1.0;

    /* Centered difference */
    for (i=1; i<nt-1; i++) {
        data[i] = (tmpdata[i+1] - tmpdata[i-1]) / (2.0 * dt);
    }

    /* Simple difference for data[0] */
    if (nt>1) data[0]    = (tmpdata[1] - tmpdata[0]) / dt;
    if (nt>2) data[nt-1] = (tmpdata[nt-1] - tmpdata[nt-2]) / dt;

    free(tmpdata);
}

/* END OF FILE */
