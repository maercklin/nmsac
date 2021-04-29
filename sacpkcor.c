/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACPKCOR - SAC cross-correlation pick refinement
*************************************************************************

Author: Nils Maercklin,
    RISSC, AMRA Scarl, Naples, Italy, December 2010

Version: 2011-06-12

Modification history:
    2010-12-16 (NM): First version of this code
    2011-02-18 (NM): Modified analysis window, pick between t1 to t2 using 
                trace samples from t1 to t2+wl; optionally save correlation 
                coefficient in SAC header (-ck option); minor other changes
    2011-02-24: Correct of pick time in case of non-zero header B
    2011-06-12: Added more flexible output file names

Notes: Input files must be evenly-sampled time series files in SAC 
    binary format (assuming NVHDR=6, IFTYPE=ITIME (1), LEVEN=1).

References:
    n/a
*************************************************************************/

#include "nmsaclib.h"     /* SAC-related functions, and other stuff */


/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACPKCOR - SAC Cross-correlation based onset time pick refinement     ",
"                                                                       ",
" Usage: sacpkcor [parameters] -f sac_files                             ",
"        sacpkcor < in_sac_file [parameters] > out_sac_file             ",
"                                                                       ",
" Parameters and defaults:                                              ",
"   -f        list of SAC binary waveform files (or read from stdin)    ",
"   -r        SAC binary file, reference waveform (or take 1st trace)   ",
"   -rk    a  SAC header field containing initial onset time pick       ",
"   -wk   t1  SAC header field for output onset time                    ",
"   -w   0.5  length of reference signal in sec (correlation template)  ",
"   -c   0.0  minimum correlation coefficient for acceptable pick       ",
"   -t1 -w/2  analysis window start (sec relative to initial pick)      ",
"   -t2 +w/2  analysis window end   (sec relative to initial pick)      ",
"                                                                       ",
"   -b1    0  Butterworth low-cut frequency in Hz  (0 = no filter)      ",
"   -b2    0  Butterworth high-cut frequency in Hz (0 = no filter)      ",
"   -bp    3  number of poles of Butterworth filters                    ",
"                                                                       ",
"   -ck none  optional SAC header field to save correlation coefficient ",
"   -a        flag: allow anti-correlation (negative corr. coefficients)",
"   -d        flag: detection, i.e. search best match on entire trace   ",
"   -o        flag: overwrite input files (default appends \".out\")    ",
"   -v        flag: verbose operation                                   ",
"                                                                       ",
"   -ext out  output file extension (if -f is set)                      ",
"   -dir      output directory (default is same as input directory)     ",
"                                                                       ",
" The picker searches for the maximum correlation coefficient between   ",
" a reference signal (template) and the current trace in a time window  ",
" from t1 to t2 (time lags relative to an initial pick). Only positive  ",
" correlation coefficients are considered, unless flag \"-a\" is set.   ",
" Initial picks, and t1, t2 are ignored in detection mode (flag -d).    ",
"                                                                       ",
" NM, 2011-06-12                                                        ",
"                                                                       ",
NULL};


/* Prototypes of functions used internally */
int get_reference_signal(SACHEAD hd, float *data, int ikey, float wl, \
    float flo, float fhi, int npoles, float **rdata);
int cc_picker(float *data, int nt, float *rdata, int rnt, int anti, float *cc);
void copy_description(SACHEAD *hd, char *rk, char *wk);

int 
main(int argc, char **argv)
{
    /* Variables */
    register int i,iarg; /* loop indices */
    SACHEAD hd;        /* SAC header */
    float *data=NULL;  /* SAC trace data */
    float *fdata=NULL; /* filtered SAC trace */
    float *rdata=NULL; /* reference (master) signal */
    float tmin;        /* start of analysis window (sec relative to rkey) */
    float tmax;        /* end of analysis window (sec relative to rkey) */
    int isfile=0;      /* file names flag */
    int isref=0;       /* internal reference trace flag */
    int is_pick=0;     /* internal valid-pick flag */
    int overwrite=0;   /* overwrite flag */
    int detect;        /* detection-mode flag */
    int anticor;       /* allow-anticorrelation flag */
    int swap;          /* internal byte-swapping flag */
    int verbose;       /* verbose flag */
    int ntr=0;         /* trace counter (for user info) */
    char *rkey=NULL;   /* header key, reference time for tmin/tmax */
    char *wkey=NULL;   /* header key, output time */
    char *ckey=NULL;   /* header key, optional output correlation coeff. */ 
    char *outfname=NULL; /* output file name */
    char *refname=NULL;/* reference trace file name */

    int irkey=0;       /* header field index of rkey */
    int iwkey=0;       /* header field index of wkey */
    int ickey=0;       /* header field index of ckey */
    float t0;          /* reference time for tmin/tmax */
    float ft, dt;      /* time of first sample and sampling rate */
    int itmin, itmax;  /* start/end of analysis window (sample indices) */
    int ntwin;         /* analysis window lengths in samples */
    int npick=0;       /* pick counter (for user info) */
    int ipick;         /* sample index of new output pick */
    float tpick=0.0;   /* time of new pick */
    float cmin=0.0;    /* min. correlation coefficient for acceptable pick */
    float c;           /* correlation coefficient */

    float wl=0.0;      /* signal time window length in seconds */
    int nwl=0;         /* time window length in samples */

    float fmin, fmax;  /* Butterworth filter frequencies */
    float f;           /* normalized frequency */
    int npoles;        /* number of poles of Butterworth filters */

    char *outdir=NULL; /* optional output directory name */
    char *outext=NULL; /* output file extension */


    /* Print documentation */
    if (argc==1) {
        i=0;
        while (sdoc[i]) fprintf(stderr, "%s\n", sdoc[i++]);
        return (EXIT_FAILURE);
    }

    /* Command line parameters */
    if (!fgetpar(argc, argv, "-w",  &wl))      wl=0.5;
    if (!fgetpar(argc, argv, "-c",  &cmin))    cmin=0.0;
    if (!fgetpar(argc, argv, "-t1", &tmin))    tmin=-0.5*wl;
    if (!fgetpar(argc, argv, "-t2", &tmax))    tmax=0.5*wl;
    if (!fgetpar(argc, argv, "-b1", &fmin))    fmin=0.0;
    if (!fgetpar(argc, argv, "-b2", &fmax))    fmax=0.0;
    if (!igetpar(argc, argv, "-bp", &npoles))  npoles=3;
    if (!sgetpar(argc, argv, "-rk", &rkey))    rkey="a";
    if (!sgetpar(argc, argv, "-wk", &wkey))    wkey="t1";
    if (!sgetpar(argc, argv, "-ck", &ckey))    ckey="NONE";

    if (!(isfile=getflag(argc, argv,  "-f")))   isfile=0;
    if (!(overwrite=getflag(argc, argv, "-o"))) overwrite=0;
    if (!(anticor=getflag(argc, argv, "-a")))   anticor=0;
    if (!(detect=getflag(argc, argv, "-d")))    detect=0;
    if (!(verbose=getflag(argc, argv, "-v")))   verbose=0;

    if (sgetpar(argc, argv, "-r", &refname))   isref=1;
    else isref=0;

    if (!sgetpar(argc, argv, "-ext", &outext)) outext="out";
    if (!sgetpar(argc, argv, "-dir", &outdir)) outdir=NULL;


    /* Get and check header field indices */
    if ((irkey=sac_hdr_index(rkey))==-1) \
        error("%s: invalid header field rk=%s\n", argv[0], rkey);
    if ((iwkey=sac_hdr_index(wkey))==-1) \
        error("%s: invalid header field wk=%s\n", argv[0], wkey);
    if (!IS_SACHDR_TMARK(iwkey)) \
        error("%s: wk=%s not allowed for output times\n", argv[0], wkey);

    ickey=sac_hdr_index(ckey);
    if (ickey>=0 && !IS_SACHDR_FLOAT(ickey)) \
        error("%s: ck=%s not allowed\n", argv[0], ckey);


    /* Check other parameters */
    cmin = fabs(cmin);
    if (wl<=0.0) \
        error("%s: window length w=%g must be positive\n", argv[0], wl);
    if (tmin>tmax) error("%s: t1 <= t2 required\n", argv[0]);
    if (!isfile && isatty(1)) \
        error("%s: can't write a SAC binary file to tty\n", argv[0]);


    /* Read reference trace */
    if (isref) {
        if (!read_sacbin_file(refname, &hd, &data, &swap)) \
            error("%s: can't read %s\n", argv[0], refname);
        if (verbose) \
            fprintf(stderr, "%s: reference file %s\n", argv[0], refname);

        nwl = get_reference_signal(hd,data,irkey,wl,fmin,fmax,npoles,&rdata);
        if (!nwl) error("%s: can't extract signal window\n", argv[0]);

        if (verbose) \
            fprintf(stderr, "%s: signal window nwl=%d samples\n", argv[0],nwl);

        free(data);
    }


    /* Loop over SAC files */
    for (iarg=isfile+1; (iarg<argc && argv[iarg][0]!='-') || !isfile; iarg++) {

        /* Read SAC binary file */
        if (isfile) {
            if (!read_sacbin_file(argv[iarg], &hd, &data, &swap)) \
                error("%s: can't read %s\n", argv[0], argv[iarg]);
            if (verbose) fprintf(stderr, "%s: file %s\n", argv[0], argv[iarg]);
        }
        else {
            if (!read_sacbin(stdin, &hd, &data, &swap)) break;
        }

        /* Read reference signal from first trace, unless read above */
        if (!isref) {
             nwl = get_reference_signal(hd, data, irkey, wl, \
                 fmin, fmax, npoles, &rdata);
             if (!nwl) error("%s: can't extract signal window\n", argv[0]);
             if (verbose) {
                 fprintf(stderr, "%s: signal window nwl=%d samples\n", \
                     argv[0], nwl);
             }
             isref = 1;
        }


        /* Get header values */
        dt = hd.delta;
        ft = hd.b;
        t0 = get_sac_hdr_float(&hd, irkey);

        /* Check required header values */
        if (ft==SAC_HEADER_FLOAT_UNDEFINED) ft = 0.0;
        if (dt<=0.0) {
            fprintf(stderr, "%s: undefined sampling rate\n", argv[0]);
            t0 = SAC_HEADER_FLOAT_UNDEFINED;     /* don't process trace */
        }

        /* Allocate space for additional, filtered trace */
        if (!(fdata=fmalloc1(hd.npts))) {
            error("%s: can't allocate space for trace, npts=%d\n", \
                argv[0], hd.npts);
        }

        /* Copy data and subtract mean */
        if (hd.depmen!=SAC_HEADER_FLOAT_UNDEFINED) {
            for (i=0; i<hd.npts; i++) fdata[i] = data[i] - hd.depmen;
        }
        else {
            memcpy(fdata, data, hd.npts*sizeof(float));
            remove_mean(fdata, hd.npts);
        }

        /* Butterworth filter(s) */
        if (npoles>0) {
            f = fmin * dt;
            if (f>0.0 && f<0.5) bfhighpass(npoles,f,hd.npts,fdata,fdata);
            f = fmax * dt;
            if (f>0.0 && f<0.5) bflowpass(npoles,f,hd.npts,fdata,fdata);
        }


        /* Analysis time window (samples), possibly cut at end of trace */
        if (detect) {
            /* entire trace... */
            itmin = 0;
            itmax = hd.npts - 1;
        }
        else {
            /* ...or window from t1 to t2+wl */
            itmin = MAX( NINT((t0+tmin-ft)/dt), 0 );
            itmax = MIN( NINT((t0+tmax+wl-ft)/dt), hd.npts-1 );
        }
        ntwin = itmax - itmin + 1;

        /* Cross-correlation analysis */
        if (ntwin>=nwl && itmin>=0 && itmax<hd.npts && \
            (t0!=SAC_HEADER_FLOAT_UNDEFINED || detect) ) {

            /* Cross-correlation picker */
            ipick = cc_picker(fdata+itmin, ntwin, rdata, nwl, anticor, &c);
            if ( (!anticor && c > cmin) || \
                 ( anticor && fabs(c)>=fabs(cmin)) ) {
                is_pick = 1;
            }
            else {
                is_pick = 0;
            }

            /* Write pick to header */
            if (is_pick) {
                npick++;

                tpick = ft + dt*(float)(itmin+ipick);
                put_sac_hdr_float(&hd, iwkey, tpick);
                copy_description(&hd, rkey, wkey);

                if (ickey>0) put_sac_hdr_float(&hd, ickey, c);
            }

            /* Diagnostic print */
            if (verbose) {
                fprintf(stderr, "%s: window [%g, %g], ", \
                    argv[0], ft+dt*(float)itmin, ft+dt*(float)itmax);
                if (is_pick) {
                    fprintf(stderr, "pick %s=%g c=%g\n", wkey, tpick, c);
                }
                else {
                    fprintf(stderr, "no pick, c=%g < cmin=%g\n", c, cmin);
                }
            }

            /* Write SAC binary file to output file... */
            if (isfile && is_pick) {
                if (overwrite) {
                    outfname=argv[iarg];
                }
                else {
                    if (!out_file_name(argv[iarg], outext, outdir, &outfname)) \
                    error("%s: can't set output file name\n", argv[0]);
                }

                if (!write_sacbin_file(outfname, hd, data, swap)) \
                    error("%s: can't write %s\n", argv[0], outfname);
                if (!overwrite) free(outfname);
            }
            /* ... or write SAC binary file to stdout */
            else if (is_pick) {
                if (!write_sacbin(stdout, hd, data, swap)) \
                    error("%s: can't write SAC trace to stdout\n", argv[0]);
            }

        }
        /* ...diagnostic print */
        else if (verbose) {
            if (t0==SAC_HEADER_FLOAT_UNDEFINED) {
                fprintf(stderr, "%s: header %s not set, trace skipped\n", \
                    argv[0], rkey);
            }
            else {
                fprintf(stderr, "%s: window [%g, %g] ", \
                    argv[0], ft+dt*(float)itmin, ft+dt*(float)itmax);
                fprintf(stderr, "too short for nwl=%d samples\n", nwl);
            }
        }

        /* Free space for trace data */
        free(data);
        free(fdata);

        /* Increment trace counter */
        ntr++;
    }

    /* Diagnostic print */
    if (verbose) {
        fprintf(stderr, "%s: picked %d of %d trace%s\n", \
            argv[0], npick, ntr, (ntr==1)?"":"s");
    }

    return EXIT_SUCCESS;
}


/************************************************************************/
/* Functions used internally                                            */
/************************************************************************/

void copy_description(SACHEAD *hd, char *rk, char *wk)
/************************************************************************
copy_description - get pick description field names and copy description
*************************************************************************
Input:
hd          SAC header  structure   (see SACHEAD typedef)
rk          header field name of time marker for reading (a,t0,...,t9)
wk          header field name of time marker for writing (a,t0,...,t9)

Output:
hd          modified SAC header structure (one of ka,kt0,...,kt9)

*************************************************************************
Author:  Nils Maercklin, 15 December 2010
*************************************************************************/
{
    char key[9];
    char pickstring[9];
    int ikey;

    sprintf(key, "k%s", rk);
    ikey = sac_hdr_index(key);

    if (ikey>=0 && get_sac_hdr_string(hd, ikey, 0, pickstring)) {
        sprintf(key, "k%s", wk);
        ikey = sac_hdr_index(key);

        if (ikey>=0) {
            put_sac_hdr_string(hd, ikey, pickstring);
        }
    }
}



int get_reference_signal(SACHEAD hd, float *data, int ikey, float wl, \
    float flo, float fhi, int npoles, float **rdata)
/************************************************************************
get_reference_signal - get reference signal for correlation picker
*************************************************************************
Input:
hd          SAC header  structure   (see SACHEAD typedef)
data        array[hd.npts] of time series
ikey        index of header field with window start time 
wl          signal time window length in seconds
flo         low-cut frequency in Hz
fhi         high-cut frequency in Hz
npoles      number of poles in Butterworth filters

Output:
rdata       array[N] of reference signal samples
            returns the number of samples N of the reference signal

Note:
Memory allocation for *rdata within this function; remember to free 
the space, if *rdata is no longer needed.

*************************************************************************
Author:  Nils Maercklin, 8 December 2010
*************************************************************************/
{
    int it, nwl, iwl;
    float dt, ft, t;
    float *fdata=NULL;

    /* Get sampling interval and window length in samples */
    dt = hd.delta;
    if (dt==SAC_HEADER_FLOAT_UNDEFINED || dt<=0.0 || wl<=0.0) {
        nwl = 0;
        return nwl;
    }
    else {
        nwl = NINT(wl/dt);
    }

    /* Index of first signal sample */
    ft = hd.b;
    t = get_sac_hdr_float(&hd, ikey);
    if (ft==SAC_HEADER_FLOAT_UNDEFINED) ft = 0.0;
    if (t==SAC_HEADER_FLOAT_UNDEFINED) {
        t = ft;
        fprintf(stderr, "%s(): %s undefined, assuming signal at t=%g\n",\
            __func__, sac_hdr_name(ikey), t);
    }
    iwl = NINT( (t - ft) / dt);

    /* Possibly adjust window length to trace length */
    if ((iwl+nwl) > hd.npts) nwl = hd.npts - iwl - 1;

    /* Process data and extract signal window */
    if (nwl>0 && iwl>=0 && iwl<hd.npts-nwl) {
        /* Allocate space for filtered trace and signal window */
        if (!(fdata=fmalloc1(hd.npts))) return 0;
        if (!(rdata[0]=fmalloc1(nwl)))  return 0;

        /* Process data */
        memcpy(fdata, data, hd.npts*sizeof(float));
        remove_mean(fdata, hd.npts);
        butterworth_bandpass(npoles, flo, fhi, dt, hd.npts, 0, fdata); 

        /* Extract window */
        for (it=0; it<nwl; it++) rdata[0][it] = fdata[it+iwl];

        /* Free space and return window length */
        free(fdata);
        return nwl;
    }
    else {
        fprintf(stderr, "%s(): window nwl=%d too long or outside range\n", \
            __func__, nwl);
        return 0;
    }
}



int cc_picker(float *data, int nt, float *rdata, int rnt, int anti, float *cc)
/************************************************************************
cc_picker - cross-correlation arrival time picker  (experimental code)
*************************************************************************
Input:
data        array[nt] of time series (seismic trace)
nt          number of samples in array data[]
rdata       array[nt] of reference signal (correlation template)
rnt         number of samples in array rdata[]
anti        flag: 1 = allow anticorrelation; 0 = allow positive CC only

Output:
cc          maximum cross-correlation coefficient
            returns the sample index (time lag) of the maximum CC

Notes:
The cross-correlation is computed from 0 to nt-rnt-1, and the length rnt 
of array rdata[] must be smaller than nt.
This is an experimental function and may be replaced by an optimized code.

*************************************************************************
Author:  Nils Maercklin, 8 December 2010
*************************************************************************/
{
    register int il, it;
    float acf, racf, ccf;
    float cmax=0.0;
    float c, rm, dm;
    int ilmax=0;


    /* Clip signal, if too long (should be avoided) */
    if (rnt>nt) rnt = nt;

    /* Mean values of time series */
    rm = dm = 0.0;
    for (it=0; it<rnt; it++) {
        rm += rdata[it];
        dm += data[it];
    }
    rm /= ((float)rnt);
    dm /= ((float)rnt);


    /* Zero-lag autocorrelation of reference signal */
    racf = 0.0;
    for (it=0; it<rnt; it++) racf += (rdata[it]-rm)*(rdata[it]-rm);
    racf = sqrt(racf);

    /* Loop over lags (from 0 to nt-rnt samples) */
    for (il=0; il<(nt-rnt); il++) {
        /* Correlation sums */
        ccf = 0.0;
        acf = 0.0;
        for (it=0; it<rnt; it++) {
            if ((il+it) < nt) {
                ccf += (data[il+it]-dm)*(rdata[it]-rm);
                acf += (data[il+it]-dm)*(data[il+it]-dm);
            }
        }
        acf = sqrt(acf);

        /* Correlation coefficient */
        if (acf!=0.0 && racf!=0.0) {
            c = ccf / (acf*racf);
        }
        else {
            c = 0.0;
        }

        /* Store maximum correlation coefficient and lag */
        if ( (anti && fabs(c)>fabs(cmax)) || (!anti && c>cmax) ) {
            cmax  = c;
            ilmax = il;
        }

        /* Update mean */
        dm -= data[il] / ((float)rnt);
        if ((il+rnt)<nt) dm += data[il+rnt] / ((float)rnt);
    }


    /* Return lag of maximum correlation coefficient */
    *cc = cmax;
    return ilmax;
}

/* END OF FILE */
