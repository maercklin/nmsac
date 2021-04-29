/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACPOFILT - SAC POlarization FILTer for three-component (3-C) data
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, February 2010

Version: 2011-06-12

Modifications:
  2010-02-04 (NM): First SAC version of this code
  2010-03-17 (NM): Added option -z for fast eigenanalysis of zero-mean 
                   data, and enabled optional zerophase bandpass
  2011-06-12 (NM): Added more flexible output file names

Notes: Input files must be evenly-sampled time series files in SAC 
    binary format (assuming NVHDR=6, IFTYPE=ITIME (1), LEVEN=1).
    Three subsequent files/traces are considered as one 3-C dataset with 
    common trace length, start time, and sampling rate.
    No station/component consistency checks are made (yet).

References:
Kanasewich, E. R. (1981). Time Sequence Analysis in Geophysics.
    The University of Alberta Press.
Montalbetti, J. R. and Kanasewich, E. R. (1970). Enhancement of 
    teleseismic body phases with a polarization filter. 
    Geophys. J. R. Astr. Soc., 21(2), 119-129.
*************************************************************************/

#include "nmsaclib.h"     /* SAC-related functions, and other stuff */


/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACPOFILT - SAC Polarization filter for three-component data          ",
"                                                                       ",
" Usage: sacpofilt [parameters] -f sac_files                            ",
"        sacpofilt <stdin [parameters] > stdout                         ",
"                                                                       ",
" Parameters and defaults:                                              ",
"   -f        list of SAC binary waveform files (or read from stdin)    ",
"   -w   0.5  correlation time window length in seconds                 ",
"   -s   0.0  smoothing time window length in seconds                   ",
"   -p    rl  rectilinearity attribute (rl, rl2, or tau)                ",
"   -q   1.0  contrast parameter of rectilinearity RL                   ",
"   -pe  1.0  exponent of rectilinearity filter function                ",
"   -de  1.0  exponent of direction filter function                     ",
"   -z        flag: assume zero mean in correlation windows (faster)    ",
"   -v        flag: verbose operation                                   ",
"                                                                       ",
"   -b1    0  Butterworth low-cut frequency in Hz  (0 = no filter)      ",
"   -b2    0  Butterworth high-cut frequency in Hz (0 = no filter)      ",
"   -bp    3  number of poles of Butterworth filters                    ",
"   -bz       flag: apply zerophase (two-pass) Butterworth filter       ",
"                                                                       ",
"   -ext      output file extension, default is \".pflt\" (if -f is set)",
"   -dir      output directory (default is same as input directory)     ",
"                                                                       ",
" This is an implementation of the Montalbetti & Kanasewich (1970)      ",
" polarization filter.                                                  ",
" Three subsequent files/traces are considered as one three-component   ",
" dataset with common trace length, start time, and sampling rate.      ",
" Default output file suffix is \".pflt\".                              ",
"                                                                       ",
" NM, 2011-06-12                                                        ",
"                                                                       ",
NULL};


/* Prototype of a function used internally */
void do_smooth(float *data, int nt, int isl);

int 
main(int argc, char **argv)
{
    /* Variables */
    register int i,j,it,iarg;/* loop indices */
    SACHEAD hd;          /* SAC header */
    SACHEAD hd3[3];      /* SAC header, three components */
    float *data=NULL;    /* SAC input trace data */
    float **data3=NULL;  /* SAC 3-C input data array */
    float *wfilt=NULL;   /* filter function, rectilinearity */
    float *dfilt=NULL;   /* filter function, direction cosines */
    float **vdata=NULL;  /* eigenvector array */
    float **ddata=NULL;  /* eigenvalue array */
    char *outfname=NULL; /* output file name */
    int isfile=0;        /* file names flag */
    int swap;            /* byte-swapping flag */
    int verbose=0;       /* verbose flag */
    int zeromean=0;      /* zero-mean flag (assume data have zero mean) */
    int nstat=0;         /* 3-C station counter (for user info) */
    int icomp=0;         /* component counter */

    int nt=0;            /* number of time samples */
    float dt;            /* time sampling interval in seconds */
    float fmin, fmax;    /* Butterworth filter frequencies */
    int npoles;          /* number of poles of Butterworth filters */
    int zerophase;       /* flag for zerophase Butterworth filters */

    float wl=0.0;        /* polarization time window length in seconds */
    float sl=0.0;        /* smoothing time window length in seconds */
    int iwl=0;           /* polarization time window length in samples */
    int isl=0;           /* smoothing time window length in samples */

    float wpow=1.0;      /* exponent of rectilinearity filter function */
    float dpow=1.0;      /* exponent of direction filter function */
    float rlq=0.0;       /* contrast factor of rectilinearity RL */
    char *pattr=NULL;    /* polarization attribute of rectilinearity */

    char *outdir=NULL;   /* optional output directory name */
    char *outext=NULL;   /* output file extension */


    /* Print documentation */
    if (argc==1) {
        i=0;
        while (sdoc[i]) fprintf(stderr, "%s\n", sdoc[i++]);
        return (EXIT_FAILURE);
    }

    /* Command line parameters */
    if (!fgetpar(argc, argv, "-w",  &wl))      wl=0.5;
    if (!fgetpar(argc, argv, "-s",  &sl))      sl=0.0;
    if (!fgetpar(argc, argv, "-q",  &rlq))     rlq=1.0;
    if (!fgetpar(argc, argv, "-pe", &wpow))    wpow=1.0;
    if (!fgetpar(argc, argv, "-de", &dpow))    dpow=1.0;
    if (!sgetpar(argc, argv, "-p",  &pattr))   pattr="rl";

    if (!fgetpar(argc, argv, "-b1", &fmin))    fmin=0.0;
    if (!fgetpar(argc, argv, "-b2", &fmax))    fmax=0.0;
    if (!igetpar(argc, argv, "-bp", &npoles))  npoles=3;
    if (!(zerophase=getflag(argc, argv, "-bz"))) zerophase=0;

    if (!(zeromean=getflag(argc, argv, "-z")) && \
        !(zeromean=getflag(argc, argv, "-zm"))) zeromean=0;

    if (!(isfile=getflag(argc, argv, "-f")))   isfile=0;
    if (!(verbose=getflag(argc, argv, "-v")))  verbose=0;

    if (!sgetpar(argc, argv, "-ext", &outext)) outext="pflt";
    if (!sgetpar(argc, argv, "-dir", &outdir)) outdir=NULL;


    /* Check parameters */
    if (fmin<=0.0 && fmax<=0.0) npoles = 0;
    if (wl<0.0) wl *= -1.0;
    if (sl<0.0) sl *= -1.0;
    if (wpow<0.0) wpow = 0.0;
    if (dpow<0.0) dpow = 0.0;
    if (!isfile && isatty(1)) \
        error("%s: can't write a SAC binary file to tty\n", argv[0]);


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


        /* Subtract mean value from trace data */
        if (hd.depmen!=SAC_HEADER_FLOAT_UNDEFINED) {
            for (it=0; it<hd.npts; it++) data[it] -= hd.depmen;
        }
        else {
            remove_mean(data, hd.npts);
        }        

        /* Butterworth low-cut and/or high-cut filter(s) */
        if (npoles>0) {
            butterworth_bandpass(npoles, fmin, fmax, dt, hd.npts, \
                zerophase, data);
        }


        /* Allocate space for three-component data */
        if (!icomp) {
            nt = hd.npts;
            if (!(data3 = fmalloc2(hd.npts, 3)) || \
                !(vdata = fmalloc2(9, hd.npts)) || \
                !(ddata = fmalloc2(3, hd.npts))) {
                error("%s: can't allocate space for 3-C data\n", argv[0]);
            }
        }

        /* Copy data */
        nt = (hd.npts>nt) ? nt : hd.npts;
        memcpy((void *)&hd3[icomp],  (const void *)&hd, sizeof(SACHEAD));
        memcpy((void *)data3[icomp], (const void *)data, nt*sizeof(float));

        /* Free space, and increment component counter */
        free(data);
        icomp++;

        /* Process three-component data */
        if (icomp==3) {
            icomp = 0;
            nstat++;

            /* Time window lengths in samples */
            iwl = NINT(wl/dt);
            if (iwl>nt-2) iwl = nt - 2;
            if (iwl<1)    iwl = 1;
            isl = NINT(sl/dt);
            if (isl>nt-2) isl = nt - 2;
            if (isl<1)    isl = 0;

            if (verbose) {
                fprintf(stderr,"%s: station %d: %d samples time window\n", \
                    argv[0], nstat, iwl);
            }


            /* Covariance analysis */
            if (zeromean) do_eigen_zm(data3, nt, iwl, ddata, vdata);
            else          do_eigen(data3, nt, iwl, ddata, vdata);


            /* Allocate space for filter traces */
            if (!(wfilt = fmalloc1(nt)) || !(dfilt = fmalloc1(nt))) {
                error("%s: can't allocate space for filter\n", argv[0]);
            }

            /* Compute rectilinearity attribute (filter trace) */
            if (wpow) {
                /* Rectilinearity RL */
                if (STRCEQ(pattr, "rl") || STRCEQ(pattr, "rl1")) {
                    for (it=0; it<nt; it++) {
                        wfilt[it] = calc_rl(ddata[it], rlq, 0);
                    }
                }
                else if (STRCEQ(pattr, "rl2")) {
                    for (it=0; it<nt; it++) {
                        wfilt[it] = calc_rl(ddata[it], rlq, 1);
                    }
                }
                else if (STRCEQ(pattr, "tau")) {
                    for (it=0; it<nt; it++) {
                        wfilt[it] = calc_tau(ddata[it]);
                    }
                }
                else {
                    fprintf(stderr, \
                        "\n%s: unsupported attribute \"%s\", ignored\n",\
                            argv[0], pattr);
                    /* Ignore attribute */
                    wpow = 0.0;
                }

                /* Smooth filter trace */
                if (isl) do_smooth(wfilt, nt, isl);

                /* Raise filter trace to power wpow */
                for (it=0; it<nt; it++) {
                    wfilt[it] = pow(wfilt[it], wpow);
                }
            }


            /* Loop over components (apply filter and write trace) */
            for (j=0; j<3; j++) {

                /* Montalbetti-Kanasewich filter */
                if (wpow) {
                    for (it=0; it<nt; it++) data3[j][it] *= wfilt[it];
                }
                if (dpow) {
                    for (it=0; it<nt; it++) dfilt[it] = fabs(vdata[it][j]);
                    if (isl) do_smooth(dfilt, nt, isl);
                    for (it=0; it<nt; it++) {
                        data3[j][it] *= pow(dfilt[it], dpow);
                    }
                }


                /* Set output SAC header */
                hd3[j].npts = nt;
                set_sac_depminmax(&hd3[j], data3[j]);


                /* Write SAC binary files... */
                if (isfile) {
                    if (!out_file_name(argv[iarg+j-2], outext, outdir, &outfname)) \
                            error("%s: can't set output file name\n", argv[0]);

                    if (!write_sacbin_file(outfname, hd3[j], data3[j], swap))\
                        error("%s: can't write %s\n", argv[0], outfname);
                    free(outfname);
                }
                /* ... or write SAC binary to stdout */
                else {
                    if (!write_sacbin(stdout, hd3[j], data3[j], swap)) \
                        error("%s: can't write SAC trace to stdout\n", argv[0]);
                }
            }

            /* Free space of three-component data and output trace */
            free2((void *) data3);
            free2((void *) vdata);
            free2((void *) ddata);
            free((void *) dfilt);
            free((void *) wfilt);
        }
    }

    /* Diagnostic print */
    if (icomp) {
        fprintf(stderr, "%s: last %d trace%s skipped\n", \
            argv[0], icomp, (icomp==1)?"":"s");
    }

    return EXIT_SUCCESS;
}


/************************************************************************/
/* Function used internally                                             */
/************************************************************************/

void do_smooth(float *data, int nt, int isl)
/************************************************************************
do_smooth - compute ellipticities (polarization attributes)

Input:
d       three-element array of eigenvalues (d[3])
i1      eigenvalue index, numerator (values 0,1,2)
i2      eigenvalue index, denominator (values 0,1,2)

Author: Nils Maercklin, 1998
*************************************************************************/
{
    register int it;
    float *tmpdata=NULL;
    float sum=0.0;


    if (isl>0 && nt>isl) {
        tmpdata=fmalloc1(nt);

        sum = 0.0;
        for (it=0; it<isl && it<nt; it++) {
            sum += data[it];
            tmpdata[it] /= ((float)(it+1));
        }
        for (it=isl/2; it<nt-isl/2; it++) {
            sum -= data[it-isl/2];
            sum += data[it+isl/2];
            tmpdata[it] = sum / ((float) isl);
        }
        for (it=nt-isl/2; it<nt; it++) {
            sum -= data[it-isl/2];
            tmpdata[it] = sum / ((float)(--isl));
        }

        memcpy((void *)data, (const void *)tmpdata, nt*sizeof(float));
        free(tmpdata);
    }
}

/* END OF FILE */
