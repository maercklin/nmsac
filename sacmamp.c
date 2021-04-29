/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACMAMP - SAC Multicomponent amplitude or RMS amplitude in time window
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, September 2009

Version: 2011-06-12

Modifications:
    2009-09-13 (NM): Initial version of this program
    2010-02-26 (NM): Added options for amplitude ratio and for energy output
    2011-06-12 (NM): Added more flexible output file names

Notes: Input files must be evenly-sampled time series files in SAC 
    binary format (assuming NVHDR=6, IFTYPE=ITIME (1), LEVEN=1).
    N subsequent files/traces are considered as one N-component dataset 
    with common trace length, start time, and sampling rate (N=1,2,3).
    No station/component consistency checks are made (yet).
*************************************************************************/

#include "nmsaclib.h"     /* SAC-related functions, and other stuff */


/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACMAMP - SAC Multicomponent Amplitude or RMS Amplitude calculation   ",
"                                                                       ",
" Usage: sacmamp [parameters] -f sac_files                              ",
"        sacmamp <stdin [parameters] > stdout                           ",
"                                                                       ",
" Parameters and defaults:                                              ",
"   -f        list of SAC binary waveform files (or read from stdin)    ",
"   -n     3  total number of components (N = 1...6)                    ",
"   -m     0  if positive, compute ratio (last N-M)/(all N) components  ",
"   -w   0.0  RMS time window length in seconds                         ",
"   -e        flag: write squared-amplitude trace (energy)              ",
"   -v        flag: verbose operation                                   ",
"                                                                       ",
"   -ext      output file extension, default is \".ampN\" (if -f is set)",
"   -dir      output directory (default is same as input directory)     ",
"                                                                       ",
" N subsequent files/traces are considered as one N-component dataset   ",
" with common trace length, start time, and sampling rate.              ",
" The default suffix for output files is \".ampN\" with N = 1...6.      ",
"                                                                       ",
" NM, 2010-02-26                                                        ",
"                                                                       ",
NULL};

#define MAXCOMP  6      /* arbitrary limit on number of components */

int 
main(int argc, char **argv)
{
    /* Variables */
    register int i,j,it,iarg;/* loop indices */
    SACHEAD hd;          /* SAC header */
    SACHEAD hd3[MAXCOMP]; /* SAC header, MAXCOMP components */
    float *data=NULL;    /* SAC trace data */
    float **data3=NULL;  /* SAC N-C input data array */
    char *outfname=NULL; /* output file name */
    char *outfext=NULL;  /* output file name extension */
    int isfile;          /* file names flag */
    int swap;            /* byte-swapping flag */
    int verbose;         /* verbose flag */
    int energy=0;        /* output-energy flag */
    int nstat=0;         /* N-C station counter (for user info) */
    int ncomp=0;         /* total number of components */
    int mcomp=0;         /* number of components in subset (for ratio) */
    int icomp=0;         /* component counter */

    float dt;            /* time sampling interval in seconds */
    int nt=0;            /* number of time samples */
    float sqrsum=0.0;    /* sum of squared amplitudes */
    float wl=0.0;        /* averaging time window length in seconds */
    int iwl=0;           /* time window length in samples */

    char *outdir=NULL;   /* optional user-specified output directory name */
    char *outext=NULL;   /* user-specified output file extension */


    /* Print documentation */
    if (argc==1) {
        i=0;
        while (sdoc[i]) fprintf(stderr, "%s\n", sdoc[i++]);
        return (EXIT_FAILURE);
    }

    /* Command line parameters */
    if (!fgetpar(argc, argv, "-w",  &wl))      wl=0.0;
    if (!igetpar(argc, argv, "-n",  &ncomp))   ncomp=3;
    if (!igetpar(argc, argv, "-m",  &mcomp))   mcomp=0;
    if (!(isfile=getflag(argc, argv, "-f")))   isfile=0;
    if (!(energy=getflag(argc, argv, "-e")))   energy=0;
    if (!(verbose=getflag(argc, argv, "-v")))  verbose=0;

    if (!sgetpar(argc, argv, "-ext", &outext)) outext=NULL;
    if (!sgetpar(argc, argv, "-dir", &outdir)) outdir=NULL;


    /* Check parameters */
    if (ncomp<1 || ncomp>MAXCOMP) \
        error("%s: -n $d: number of components must be in [1, %d]\n", \
            argv[0], ncomp, MAXCOMP);
    if (mcomp<0 || mcomp>=ncomp) \
        error("%s: -m $d: number of components must be in [0, %d]\n", \
            argv[0], mcomp, ncomp-1);

    if (!isfile && isatty(1)) \
        error("%s: can't write a SAC binary file to tty\n", argv[0]);


    /* Set user-secified output file extension */
    if (outext) {
        outfext=outext;
    }
    /* ...or set default extension */
    else {
        outfext=malloc(5);
        if (mcomp) sprintf(outfext, "r%dt%d", ncomp-mcomp, ncomp);
        else sprintf(outfext, "amp%d", ncomp);
    }


    /* Loop over SAC files */
    icomp = 0;
    nstat = 0;
    for (iarg=isfile+1; (iarg<argc && argv[iarg][0]!='-') || !isfile; iarg++) {

        /* Read SAC binary file */
        if (isfile) {
            if (!read_sacbin_file(argv[iarg], &hd, &data, &swap)) \
                error("%s: can't read %s\n", argv[0], argv[iarg]);
            if (verbose) \
                fprintf(stderr, "%s: station %d, file %d: %s\n", \
                    argv[0], nstat+1, icomp+1, argv[iarg]);
        }
        else {
            if (!read_sacbin(stdin, &hd, &data, &swap)) break;
        }

        /* Get required header values */
        dt = hd.delta;

        /* Check header values */
        if (dt<=0.0) {
            fprintf(stderr,"%s: invalid delta=%g, assuming delta=1.0\n", \
                argv[0], hd.delta);
            dt = 1.0;
        }


        /* Allocate space for N-component data */
        if (!icomp) {
            nt = hd.npts;
            if (!( data3 = fmalloc2(hd.npts, ncomp))) {
                error("%s: can't allocate space for %d-C data\n", \
                    argv[0], ncomp);
            }
        }

        /* Copy data */
        nt = (hd.npts>nt) ? nt : hd.npts;
        memcpy((void *)&hd3[icomp], (const void *)&hd, sizeof(SACHEAD));
        memcpy((void *)data3[icomp], (const void *)data, nt*sizeof(float));

        /* Free space, and increment component counter */
        free(data);
        icomp++;

        /* Process three-component data */
        if (icomp==ncomp) {
            icomp = 0;
            nstat++;

            /* Allocate space for output data */
            if (!(data=fmalloc1(nt))) \
                error("%s: can't allocate space for output trace\n", argv[0]);

            /* Get time window length */
            iwl = (wl>dt) ? NINT(wl/dt) : 1;

            if (verbose && iwl>1) {
                fprintf(stderr,"%s: station %d: %d samples time window\n", \
                    argv[0], nstat, iwl);
            }

            /* Initialize running sum */
            sqrsum = 0.0;

            /* Loop over samples, total multicomponent energy */
            for (it=0; it<nt; it++) {
                /* Time window */
                if (it >= iwl) {
                    for (j=0; j<ncomp; j++) \
                        sqrsum -= data3[j][it-iwl]*data3[j][it-iwl];
                }

                /* Running average in time window */
                for (j=0; j<ncomp; j++) sqrsum += data3[j][it]*data3[j][it];

                /* Total energy trace */
                data[it] = sqrsum / ((float)iwl);
            }

            /* Loop over samples (energy ratio) */
            if (mcomp) {
                sqrsum = 0.0;

                for (it=0; it<nt; it++) {
                    /* Time window */
                    if (it >= iwl) {
                        for (j=mcomp; j<ncomp; j++) \
                            sqrsum -= data3[j][it-iwl]*data3[j][it-iwl];
                    }

                    /* Running average in time window */
                    for (j=mcomp; j<ncomp; j++) \
                        sqrsum += data3[j][it]*data3[j][it];

                    /* Energy ratio trace */
                    if (data[it]) data[it] = sqrsum / (data[it] * (float)iwl);
                    else          data[it] = 0.0;

                }
            }

            /* Output amplitude (ratio) trace */
            if (!energy) {
                for (it=0; it<nt; it++) data[it] = sqrt(data[it]);
            }


            /* Set output SAC header */
            hd3[0].npts = nt;
            set_sac_depminmax(&hd3[0], data);
            put_sac_hdr_string(&hd3[0], sac_hdr_index("kcmpnm"), outfext);


            /* Write SAC binary files... */
            if (isfile) {
                if (!out_file_name(argv[iarg-ncomp+1], outfext, outdir, &outfname)) \
                    error("%s: can't set output file name\n", argv[0]);

                if (!write_sacbin_file(outfname, hd3[0], data, swap))\
                    error("%s: can't write %s\n", argv[0], outfname);
                free(outfname);
            }
            /* ... or write SAC binary to stdout */
            else {
                if (!write_sacbin(stdout, hd3[0], data, swap)) \
                    error("%s: can't write SAC trace to stdout\n", argv[0]);
            }

            /* Free space of output and three-component data */
            free2((void *) data3);
            free((void *) data);
        }
    }

    /* Diagnostic print */
    if (icomp) {
        fprintf(stderr, "%s: last %d trace%s skipped\n", \
            argv[0], icomp, (icomp==1)?"":"s");
    }

    return EXIT_SUCCESS;
}

/* END OF FILE */
