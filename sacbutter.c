/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACBUTTER - SAC Butterworth high-pass and/or low-pass filter
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, September 2009

Version: 2011-06-12

Modifications:
    2009-09-13 (NM): First version of this code 
    2011-06-12 (NM): Added more flexible output file names 

Notes: Input files must be evenly-sampled time series files in SAC 
    binary format (assuming NVHDR=6, IFTYPE=ITIME (1), LEVEN=1).
*************************************************************************/

#include "nmsaclib.h"     /* SAC-related functions, and other stuff */

/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACBUTTER - SAC Butterworth high-pass and/or low-pass filter          ",
"                                                                       ",
" Usage: sacbutter [parameters] -f sac_files                            ",
"        sacbutter <stdin [parameters] > stdout                         ",
"                                                                       ",
" Parameters and defaults:                                              ",
"   -f        list of SAC binary waveform files (or read from stdin)    ",
"                                                                       ",
"   -b1    0  Butterworth low-cut frequency in Hz  (0 = no filter)      ",
"   -b2    0  Butterworth high-cut frequency in Hz (0 = no filter)      ",
"   -bp    3  number of poles of Butterworth filters                    ",
"                                                                       ",
"   -bz       flag: apply zerophase (two-pass) filter                  ",
"   -o        flag: overwrite input files (default appends \".flt\")    ",
"                                                                       ",
"   -ext flt  output file extension (if -f is set)                      ",
"   -dir      output directory (default is same as input directory)     ",
"                                                                       ",
" The parameters -flo and -fhi are acceptable aliases for -b1 and -b2,  ", 
" respectively. If the overwrite flag -o is set, original data are lost.",
"                                                                       ",
" NM, 2011-06-12                                                        ",
"                                                                       ",
NULL};

/* Definitions */
#define RAD2DEG  57.29577951

int 
main(int argc, char **argv)
{
    /* Variables */
    register int iarg, it; /* loop indices (files, time samples) */
    SACHEAD hd;          /* SAC header */
    float *data=NULL;    /* SAC input trace data */
    char *outfname=NULL; /* output file name */
    int isfile=0;        /* file names flag */
    int isvalid=0;       /* internal SAC trace validation flag */
    int swap;            /* byte-swapping flag */
    int overwrite=0;     /* overwrite flag */
    int zerophase=0;     /* zero-phase flag */

    float fmin, fmax;    /* Butterworth filter frequencies */
    int npoles;          /* number of poles of Butterworth filters */

    char *outdir=NULL;   /* optional output directory name */
    char *outext=NULL;   /* output file extension */


    /* Print documentation */
    if (argc==1) {
        int i=0;
        while (sdoc[i]) fprintf(stderr, "%s\n", sdoc[i++]);
        return (EXIT_FAILURE);
    }

    /* Command line parameters */
    if (!fgetpar(argc, argv, "-flo", &fmin) && \
        !fgetpar(argc, argv, "-b1",  &fmin))      fmin=0.0;
    if (!fgetpar(argc, argv, "-fhi", &fmax) && \
        !fgetpar(argc, argv, "-b2",  &fmax))      fmax=0.0;
    if (!igetpar(argc, argv, "-bp",  &npoles))    npoles=3;

    if (!(zerophase=getflag(argc, argv, "-z")) && \
        !(zerophase=getflag(argc, argv, "-bz")))  zerophase=0;
    if (!(overwrite=getflag(argc, argv, "-o")))   overwrite=0;
    if (!(isfile=getflag(argc, argv, "-f")))      isfile=0;

    if (!sgetpar(argc, argv, "-ext", &outext))    outext="flt";
    if (!sgetpar(argc, argv, "-dir", &outdir))    outdir=NULL;


    /* Check parameters */
    if (npoles<=0) error("%s: number of poles must be positive\n", argv[0]);
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

        /* Check for valid trace type and sampling rate */
        isvalid = (!hd.leven || hd.iftype==IRLIM || hd.iftype==IAMPH) ? 0 : 1;
        if (hd.delta<=0.0) {
            fprintf(stderr,"%s: invalid delta=%g, skipping %s\n", \
                argv[0], hd.delta, (isfile) ? argv[iarg] : "stdin");
            isvalid = 0;
        }

        if (isvalid) {
            /* Subtract mean value from trace */
            if (hd.depmen!=SAC_HEADER_FLOAT_UNDEFINED) {
                for (it=0; it<hd.npts; it++) data[it] -= hd.depmen;
            }
            else {
                remove_mean(data, hd.npts);
            }

            /* Butterworth low-cut and/or high-cut filter(s) */
            if (npoles>0) {
                butterworth_bandpass(npoles, fmin, fmax, hd.delta, hd.npts, \
                    zerophase, data);
            }


            /* Re-compute depmen/depmin/depmax */
            set_sac_depminmax(&hd, data);

            /* Write SAC binary files... */
            if (isfile) {
                if (overwrite) {
                    outfname = argv[iarg];
                }
                else {
                    if (!out_file_name(argv[iarg], outext, outdir, &outfname)) \
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
        }
    }

    return EXIT_SUCCESS;
}

/* END OF FILE */
