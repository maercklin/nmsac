/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACROT3 - SAC Rotation of three-component data
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, February 2010

Version: 2011-06-12

Modifications:
    2010-01-27 (NM): Initial version of this program
    2010-02-26 (NM): Fixed signs in rotation ZNE -> LQT
    2011-06-12 (NM): Added more flexible output file names 

Notes: Input files must be evenly-sampled time series files in SAC 
    binary format (assuming NVHDR=6, IFTYPE=ITIME (1), LEVEN=1).
    Three subsequent files/traces are considered as one 3-C dataset with 
    common trace length, start time, and sampling rate (order: ZNE).
    No station/component consistency checks are made (yet).
    The SAC header fields CMPAZ and CMPINC are adjusted, but not KCMPNM.

Reference:
    n/a
*************************************************************************/

#include "nmsaclib.h"     /* SAC-related functions, and other stuff */


/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACROT3 - SAC Rotation of three-component data                        ",
"                                                                       ",
" Usage: sacrot3 [parameters] -f sac_files                              ",
"        sacrot3 <stdin [parameters] > stdout                           ",
"                                                                       ",
" Parameters and defaults:                                              ",
"   -f        list of SAC binary waveform files (or read from stdin)    ",
"   -a   BAZ  horizontal rotation angle in degrees (e.g. backazimuth)   ",
"   -i   0.0  vertical rotation angle in degrees (e.g. incidence angle) ",
"   -h        flag: horizontal rotation of two-component data           ",
"   -o        flag: overwrite input files (default appends \".rot\")    ",
"   -v        flag: verbose operation                                   ",
"                                                                       ",
"   -ext rot  output file extension (if -f is set)                      ",
"   -dir      output directory (default is same as input directory)     ",
"                                                                       ",
" Three subsequent files/traces are considered as one three-component   ",
" dataset with common trace length, start time, and sampling rate.      ",
" If option \"-h\" is set, only two subsequent files are considered as  ",
" one dataset of two horizontal components.                             ",
" The first component of each three-component dataset is assumed to     ",
" be the vertical component. Default output file suffix is \".rot\".    ",
"                                                                       ",
" NM, 2011-06-12                                                        ",
"                                                                       ",
NULL};


int 
main(int argc, char **argv)
{
    /* Variables */
    register int i,j,it,iarg;/* loop indices */
    SACHEAD hd;         /* SAC header, current input trace */
    SACHEAD hd3[3];     /* SAC header, three components */
    float *data=NULL;   /* SAC trace data */
    float **data3=NULL; /* SAC 3-C input data array */
    float **odata3=NULL;/* SAC 3-C output data array */
    char *outfname=NULL;/* output file name */
    int isfile;         /* file names flag */
    int swap;           /* byte-swapping flag */
    int verbose=0;      /* verbose flag */
    int overwrite=0;    /* overwrite flag */
    int geom=0;         /* geometry flag (for horizontal rotation, BAZ) */
    int horiz=0;        /* flag for horizontal rotation of 2-C data */
    int ncomp=0;        /* number of components to read (2 or 3) */
    int icomp=0;        /* component counter */
    int nstat=0;        /* 3-C station counter (for user info) */
    int nt=0;           /* number of time samples */
    float inc=0.0;      /* vertical rotation angle (incidence angle) */
    float azi=0.0;      /* horizontal rotation angle (e.g. backazimuth) */
    float theta, beta;  /* rotation angles in radians */

    char *outdir=NULL;   /* optional output directory name */
    char *outext=NULL;   /* output file extension */


    /* Print documentation */
    if (argc==1) {
        i=0;
        while (sdoc[i]) fprintf(stderr, "%s\n", sdoc[i++]);
        return (EXIT_FAILURE);
    }

    /* Command line parameters */
    geom=0; azi=0.0;
    if (!fgetpar(argc, argv, "-a", &azi))        geom=1;
    if (!fgetpar(argc, argv, "-i", &inc))        inc=0.0;

    if (!(isfile=getflag(argc, argv, "-f")))     isfile=0;
    if (!(horiz=getflag(argc, argv, "-h")))      horiz=0;
    if (!(overwrite=getflag(argc, argv, "-o")))  overwrite=0;
    if (!(verbose=getflag(argc, argv, "-v")))    verbose=0;

    if (!sgetpar(argc, argv, "-ext", &outext))   outext="rot";
    if (!sgetpar(argc, argv, "-dir", &outdir))   outdir=NULL;


    /* Check parameters */
    if (!isfile && isatty(1)) \
        error("%s: can't write a SAC binary file to tty\n", argv[0]);

    /* Number of components to read (2 or 3) */
    ncomp = (horiz) ? 2 : 3;

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


        /* Allocate space for two- or three-component data */
        if (!icomp) {
            nt = hd.npts;
            if (!( data3 = fmalloc2(hd.npts, ncomp)) || \
                !(odata3 = fmalloc2(hd.npts, ncomp))) {
                error("%s: can't allocate space for %d-C data\n", \
                    ncomp, argv[0]);
            }
        }

        /* Copy data */
        nt = (hd.npts>nt) ? nt : hd.npts;
        memcpy((void *)&hd3[icomp], (const void *)&hd, sizeof(SACHEAD));
        memcpy((void *)data3[icomp], (const void *)data, nt*sizeof(float));

        /* Free space, and increment component counter */
        free(data);
        icomp++;

        /* Process two- or three-component data */
        if (icomp==ncomp) {
            icomp = 0;
            nstat++;

            /* Get angle azi from header */
            if (geom) {
                if (hd3[ncomp-2].baz  ==SAC_HEADER_FLOAT_UNDEFINED || \
                    hd3[ncomp-2].cmpaz==SAC_HEADER_FLOAT_UNDEFINED) {
                    fprintf(stderr, "%s: station %d: BAZ or CMPAZ undefined\n",\
                        argv[0], nstat);
                    azi = 0.0;
                }
                else {
                    azi = hd3[ncomp-2].baz - hd3[ncomp-2].cmpaz;
                }
            }

            /* Rotation angles in radians */
            theta = inc * PI/180.0;
            beta  = azi * PI/180.0;

            /* Echo rotation angles */
            if (verbose) {
                if (horiz) fprintf(stderr, "%s: station %d: a=%g deg\n", \
                    argv[0], nstat, azi);
                else fprintf(stderr, "%s: station %d: a=%g deg i=%g deg\n", \
                    argv[0], nstat, azi, inc);
            }


            /* Loop over samples (2-C horizontal rotation) */
            if (horiz) {
                for (it=0; it<nt; it++) {
                    /* NE -> R */
                    odata3[0][it] =   data3[0][it] * cos(beta) \
                                    + data3[1][it] * sin(beta);
                    /* NE -> T */
                    odata3[1][it] = - data3[0][it] * sin(beta) \
                                    + data3[1][it] * cos(beta);
                }
            }
            /* Loop over samples (general 3-C rotation) */
            else {
                for (it=0; it<nt; it++) {
                    /* ZNE -> L or Z */
                    odata3[0][it] = data3[0][it] * cos(theta) - \
                                    data3[1][it] * sin(theta)*cos(beta) - \
                                    data3[2][it] * sin(theta)*sin(beta);
                    /* ZNE -> Q or R */
                    odata3[1][it] = data3[0][it] * sin(theta) + \
                                    data3[1][it] * cos(theta)*cos(beta) + \
                                    data3[2][it] * cos(theta)*sin(beta);
                    /* ZNE -> Q or R */
                    odata3[2][it] = - data3[1][it] * sin(beta) + \
                                      data3[2][it] * cos(beta);
                }
            }

            /* Output */
            for (j=0; j<ncomp; j++) {
                /* Set header */
                hd3[j].npts = nt;
                if (!horiz && j<2 && \
                    hd3[j].cmpinc!=SAC_HEADER_FLOAT_UNDEFINED) {
                    hd3[j].cmpinc += inc;
                }

                if (hd3[j].cmpaz!=SAC_HEADER_FLOAT_UNDEFINED) {
                    hd3[j].cmpaz = fmod((hd3[j].cmpaz + azi), 360.0);
                    if (hd3[j].cmpaz < 0.0) hd3[j].cmpaz += 360.0;
                }

                /* Re-compute depmen/depmin/depmax */
                if (horiz) {
                    if (azi!=0.0) set_sac_depminmax(&hd3[j], odata3[j]);
                }
                else {
                    if (j==0 && inc!=0.0) 
                        set_sac_depminmax(&hd3[j], odata3[j]);
                    else if (j && azi!=0.0) 
                        set_sac_depminmax(&hd3[j], odata3[j]);
                }


                /* Write SAC binary files... */
                if (isfile) {
                    if (overwrite) {
                        outfname = argv[iarg+j-ncomp+1];
                    }
                    else {
                        if (!out_file_name(argv[iarg+j-ncomp+1], outext, outdir, &outfname)) \
                            error("%s: can't set output file name\n", argv[0]);
                    }

                    if (!write_sacbin_file(outfname, hd3[j], odata3[j], swap))\
                        error("%s: can't write %s\n", argv[0], outfname);
                    if (!overwrite) free(outfname);
                }
                /* ... or write SAC binary to stdout */
                else {
                    if (!write_sacbin(stdout, hd3[j], odata3[j], swap)) \
                        error("%s: can't write SAC trace to stdout\n", argv[0]);
                }
            }

            /* Free space of three-component data */
            free2((void *) data3);
            free2((void *)odata3);
        }
    }

    /* Diagnostic print */
    if (icomp) {
        fprintf(stderr, "%s: last %d trace%s skipped\n", \
            argv[0], icomp, (icomp==1)?"":"s");
    }
    if (verbose) {
        fprintf(stderr, "%s: processed %d %d-C seismogram%s\n", \
            argv[0], nstat, ncomp, (nstat==1)?"":"s");
    }

    return EXIT_SUCCESS;
}

/* END OF FILE */
