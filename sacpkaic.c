/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACPKAIC - SAC single-trace onset time estimation based on AIC
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, September 2008

Version: 2011-07-04

Modification history:
    2008-03-04 (NM): Initial coding for SEG-Y traces and CWP/SU libs
    2008-09-10 (NM): First version for SAC binary traces using "nmsaclib.c",
                modified filters
    2008-10-06 (NM): Read multiple traces from stdin
    2010-06-02 (NM): Changed analysis window (picking only from t1 to t2, AR
                starts at t1-tc), added minimum SNR for acceptable pick, 
                and added experimental pick quality string
    2011-06-12 (NM): Added more flexible output file names 
    2011-07-04 (NM): Added optional output of pick error (experimental) 

Notes: Input files must be evenly-sampled time series files in SAC 
    binary format (assuming NVHDR=6, IFTYPE=ITIME (1), LEVEN=1).

References:
T. Kvaerna (1994). Accurate determination of phase arrival times 
    using autoregressive likelihood estimation. Annali di Geofisica,
    37(3), 287-300, http://hdl.handle.net/2122/1858.
M. Leonard and B.L.N. Kennett (1999). Multi-component autoregressive
    techniques for the analysis of seismograms. Phys. Earth Planet.
    Int., 113, 247-263, doi:10.1016/S0031-9201(99)00054-0.
R. Sleeman and T. van Eck (1999). Robust automatic P-phase picking:
    an on-line implementation in the analysis of broadband seismogram
    recordings. Phys. Earth Planet. Int., 113, 265-275, 
    doi:10.1016/S0031-9201(99)00007-2.
*************************************************************************/

#include "nmsaclib.h"     /* SAC-related functions, and other stuff */


/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACPKAIC - SAC single-trace onset time estimation based on AIC        ",
"            incl. optional autoregressive and/or frequency filtering   ",
"                                                                       ",
" Usage: sacpkaic [parameters] -f sac_files                             ",
"        sacpkaic < in_sac_file [parameters] > out_sac_file             ",
"                                                                       ",
" Parameters and defaults:                                              ",
"   -f        list of SAC binary waveform files (or read from stdin)    ",
"   -rk    a  SAC header field containing initial onset time pick       ",
"   -wk   t1  SAC header field for output onset time                    ",
"   -t1 -2.0  analysis window start (sec relative to initial pick)      ",
"   -t2  2.0  analysis window end   (sec relative to initial pick)      ",
"                                                                       ",
"   -tc  2.0  window length (sec) for determination of AR coefficients  ",
"   -m     0  number of AR filter coefficients (0 = no filter)          ",
"                                                                       ",
"   -b1    0  Butterworth low-cut frequency in Hz  (0 = no filter)      ",
"   -b2    0  Butterworth high-cut frequency in Hz (0 = no filter)      ",
"   -bp    3  number of poles of Butterworth filters                    ",
"                                                                       ",
"   -snr 1.0  minimum signal-to-noise ratio for acceptable pick         ",
"                                                                       ",
"   -qp    ?  phase letter in pick description field (pick quality)     ",
"   -qf .025  fraction of AIC range to estimate error (pick quality)    ",
"   -ek none  optional SAC header field to save estimated pick error    ",
"                                                                       ",
"   -x     o  output option: o=original data, r=raw data in window,     ",
"                            f=filtered data in window, a=AIC trace     ",
"   -o        flag: overwrite input files (default appends \".out\")    ",
"   -v        flag: verbose operation                                   ",
"                                                                       ",
"   -ext out  output file extension (if -f is set)                      ",
"   -dir      output directory (default is same as input directory)     ",
"                                                                       ",
" Picking is done in the time window from t1 to t2, and the optional    ",
" AR filter coefficients are determined in a window from t1-tc to t1.   ",
" Pick quality is only evaluated, if the parameter -qp is set.          ",
"                                                                       ",
" NM, 2011-06-12                                                        ",
"                                                                       ",
NULL};


/* Prototypes of functions used internally */
void pefilter(float *data, int nt, int istart, int ncorr, int nlag, \
    float *pedata);
int compute_aick(float *ndata, float *sdata, int n, int istart, \
    float *snr, float *aic);
int aic_pick_quality(float *aic, int n, float dt, float snr, float frac, \
    float wscale, char *phase, char *pickstring, float *terr);

int 
main(int argc, char **argv)
{
    /* Variables */
    register int i,iarg; /* loop index */
    SACHEAD hd;        /* SAC header */
    float *data=NULL;  /* SAC trace data */
    float tmin;        /* start of analysis window (sec relative to rkey) */
    float tmax;        /* end of analysis window (sec relative to rkey) */
    float tcorr;       /* window length for AR coeff. determination (sec) */
    int mop;           /* number of AR coefficients, i.e. operator length */
    int isfile;        /* file names flag */
    int overwrite;     /* overwrite flag */
    int swap;          /* byte-swapping flag */
    int verbose;       /* verbose flag */
    int ntr=0;         /* trace counter (for user info) */
    char *rkey=NULL;   /* header key, reference time for tmin/tmax */
    char *wkey=NULL;   /* header key, output time */
    char *ekey=NULL;   /* header key, optional output pick error */ 
    char *outfname=NULL; /* output file name */
    char *outopt=NULL; /* output option (original or windowed trace) */
    float *wdata=NULL; /* data in analysis time window (incl. tcorr win.) */
    float *fdata=NULL; /* filtered data in analysis window */
    float *adata=NULL; /* AIC function in analysis window */
    float *odata=NULL; /* pointer to output data array (e.g. data) */

    int irkey=0;       /* header field index of rkey */
    int iwkey=0;       /* header field index of wkey */
    int iekey=0;       /* header field index of ekey */
    float t0;          /* reference time for tmin/tmax */
    float ft, dt;      /* time of first sample and sampling rate */
    int itmin, itmax;  /* start/end of analysis window (sample indices) */
    int ncorr, ntwin;  /* window lengths in samples */
    int npick=0;       /* pick counter (for user info) */
    int ipick;         /* sample index of new output pick */
    float tpick=0.0;   /* time of new pick */
    float snr=0.0;     /* signal-to-noise ratio */
    float minsnr=0.0;  /* minimum signal-to-noise ratio for acceptable pick */

    float fmin, fmax;  /* Butterworth filter frequencies */
    float f;           /* normalized frequency */
    int npoles;        /* number of poles of Butterworth filters */

    char *outdir=NULL; /* optional output directory name */
    char *outext=NULL; /* output file extension */

    char pickstring[8];/* pick quality string for description field */
    char qkey[8];      /* header field name for quality string */
    char *qphase=NULL; /* phase name in quality string, e.g. P or S */
    float qfrac=0.0;   /* fraction of AIC range for quality control */
    float qscale=1.0;  /* scale factor for error->weight table */ 
    float terr=0.0;    /* estimated picking error in seconds */
    int iqkey=-1;      /* header field index of qkey */
    int quality=0;     /* quality evaluation flag */



    /* Print documentation */
    if (argc==1) {
        i=0;
        while (sdoc[i]) fprintf(stderr, "%s\n", sdoc[i++]);
        return (EXIT_FAILURE);
    }

    /* Command line parameters */
    if (!fgetpar(argc, argv, "-t1", &tmin))    tmin=-2.0;
    if (!fgetpar(argc, argv, "-t2", &tmax))    tmax=2.0;
    if (!fgetpar(argc, argv, "-tc", &tcorr))   tcorr=2.0;
    if (!fgetpar(argc, argv, "-snr",&minsnr))  minsnr=1.0;
    if (!fgetpar(argc, argv, "-b1", &fmin))    fmin=0.0;
    if (!fgetpar(argc, argv, "-b2", &fmax))    fmax=0.0;
    if (!igetpar(argc, argv, "-bp", &npoles))  npoles=3;
    if (!igetpar(argc, argv, "-m",  &mop))     mop=0;
    if (!sgetpar(argc, argv, "-rk", &rkey))    rkey="a";
    if (!sgetpar(argc, argv, "-wk", &wkey))    wkey="t1";
    if (!sgetpar(argc, argv, "-ek", &ekey))    ekey="NONE";
    if (!sgetpar(argc, argv, "-x",  &outopt))  outopt="o";
    if (!sgetpar(argc, argv, "-qp", &qphase))  qphase="?";
    if (!fgetpar(argc, argv, "-qf", &qfrac))   qfrac=0.025;
    if (!fgetpar(argc, argv, "-qs", &qscale))  qscale=1.0;

    if (!(isfile=getflag(argc, argv,  "-f")))   isfile=0;
    if (!(overwrite=getflag(argc, argv, "-o"))) overwrite=0;
    if (!(quality=getflag(argc, argv, "-qp")))  quality=0;
    if (!(verbose=getflag(argc, argv, "-v")))   verbose=0;

    if (!sgetpar(argc, argv, "-ext", &outext)) outext="out";
    if (!sgetpar(argc, argv, "-dir", &outdir)) outdir=NULL;


    /* Get and check header field indices */
    if ((irkey=sac_hdr_index(rkey))==-1) \
        error("%s: invalid header field rk=%s\n", argv[0], rkey);
    if ((iwkey=sac_hdr_index(wkey))==-1) \
        error("%s: invalid header field wk=%s\n", argv[0], wkey);
    if (!IS_SACHDR_TMARK(iwkey)) \
        error("%s: wk=%s not allowed\n", argv[0], wkey);

    iekey=sac_hdr_index(ekey);
    if (iekey>=0) {
        quality = 1;
        if (!IS_SACHDR_FLOAT(iekey)) \
            error("%s: ek=%s not allowed\n", argv[0], ekey);
    }
    if (quality) {
        sprintf(qkey, "k%s", wkey);
        if ((iqkey=sac_hdr_index(qkey))==-1) quality = 0;
    }

    /* Check other parameters */
    if (tmin>=tmax) error("%s: t1 < t2 required\n", argv[0]);
    if (mop && tcorr<=0.0) error("%s: AR filter requires tc > 0\n", argv[0]);
    if (!isfile && isatty(1)) \
        error("%s: can't write a SAC binary file to tty\n", argv[0]);
    if (mop<0) mop = abs(mop);
    if (!mop) tcorr= 0.0;

    /* Protect original data in case of windowed output */
    if (overwrite && outopt[0]!='o') overwrite=0;


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

        /* Analysis time window (samples), possibly cut at end of trace */
        itmin = MAX( NINT((t0+tmin-tcorr-ft)/dt), 0 );
        itmax = MIN( NINT((t0+tmax-ft)/dt), hd.npts-1 );
        ntwin = itmax - itmin + 1;

        /* Correlation window length for AR coefficients */
        if (mop) {
            ncorr = NINT(tcorr/dt);
            if (ncorr<=mop) ncorr = mop + 1; /* ensure appropriate length */
        }
        else {
            ncorr = 0;
        }

        /* Processing, if window inside trace etc... */
        if (ntwin>=ncorr && itmin>=0 && itmax<hd.npts && \
            t0!=SAC_HEADER_FLOAT_UNDEFINED) {

            /* Allocate space */
            wdata = fmalloc1(ntwin);
            adata = fmalloc1(ntwin);

            /* Extract time window */
            memcpy(wdata, data+itmin, ntwin*sizeof(float));

            /* Subtract mean in time window */
            remove_mean(wdata, ntwin);

            /* Autoregressive filter */
            if (mop && mop<ncorr) {
                fdata = fmalloc1(ntwin);
                pefilter(wdata, ntwin, 0, ncorr, mop, fdata);
            }
            else {
                fdata = wdata;
            }

            /* Butterworth filter(s) */
            if (npoles>0) {
                f = fmin * dt;
                if (f>0.0 && f<0.5) bfhighpass(npoles,f,ntwin,fdata,fdata);
                f = fmax * dt;
                if (f>0.0 && f<0.5) bflowpass(npoles,f,ntwin,fdata,fdata);
            }

            /* AIC and picking */
            ipick = compute_aick(fdata, fdata, ntwin, ncorr, &snr, adata);
            if (snr<minsnr) ipick = 0;
            tpick = (ipick) ? ft + dt*((float)(itmin+ipick)) : t0;

            /* Pick quality (experimental) */
            if (ipick && quality) {
                if (aic_pick_quality(adata, ntwin, dt, snr, qfrac, qscale, \
                        qphase, pickstring, &terr)) {
                    put_sac_hdr_string(&hd, iqkey, pickstring);

                    if (iekey>0) put_sac_hdr_float(&hd, iekey, terr);
                }
                else {
                    ipick = 0;
                }
            }

            /* Diagnostic print */
            if (verbose) {
                fprintf(stderr, "%s: window [%g, %g], ", \
                    argv[0], ft+dt*(float)itmin, ft+dt*(float)itmax);
                if (ipick) {
                    fprintf(stderr, "pick %s=%g snr=%g\n", wkey, tpick, snr);
                }
                else {
                    fprintf(stderr, "no pick, snr=%g\n", snr);
                }
            }

            /* Write pick time to header, increment pick counter */
            if (ipick) {
                put_sac_hdr_float(&hd, iwkey, tpick);
                npick++;
            }

            /* Select trace data for (windowed) output */
            switch (outopt[0]) {
                case 'r':
                    odata=data+itmin;
                    break;
                case 'f':
                    odata=fdata;
                    break;
                case 'a':
                    odata=adata;
                    break;
                case 'o':
                default:
                    odata=data;
                    break;
            }

            /* Set header, if windowed output requested */
            if (outopt[0]=='r' || outopt[0]=='f' || outopt[0]=='a') {
                hd.npts = ntwin;
                hd.b    = ft + dt*((float)(itmin));
                hd.depmen = 0.0;
                set_sac_depminmax(&hd, odata);
            }

            /* Write SAC binary file to output file... */
            if (isfile) {
                if (overwrite) {
                    outfname=argv[iarg];
                }
                else {
                    if (!out_file_name(argv[iarg], outext, outdir, &outfname)) \
                    error("%s: can't set output file name\n", argv[0]);
                }

                if (!write_sacbin_file(outfname, hd, odata, swap)) \
                    error("%s: can't write %s\n", argv[0], outfname);
                if (!overwrite) free(outfname);
            }
            /* ... or write SAC binary file to stdout */
            else {
                if (!write_sacbin(stdout, hd, odata, swap)) \
                    error("%s: can't write SAC trace to stdout\n", argv[0]);
            }

            /* Free space for time window and AIC */
            free(wdata);
            free(adata);
            if (mop && mop<ncorr) free(fdata);
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
                fprintf(stderr, "too short for ncorr=%d samples\n", ncorr);
            }
        }

        /* Free space for trace data */
        free(data);

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

int compute_aick(float *ndata, float *sdata, int n, int istart, \
    float *snr,float *aic)
/************************************************************************
compute_aick - compute Akaike Information Criterion (AIC) for separation
               point K between two stationary time sequences
*************************************************************************
Input:
ndata       array[n] of time series representing "noise"
sdata       array[n] of time series representing "signal"
n           number of samples in ndata[] and sdata[]
istart      index of first sample for AIC computation

Output:
            returns sample index of AIC minimum (0 = failure)
snr         signal-to-noise ratio
aic         array[n] of AIC values

*************************************************************************
References:
H. Akaike (1973). Information theory and an extension of the maximum 
    likelihood principle. In: B.N. Petrov and F. Caski (editors), 
    Proceedings of the Second International Symposium on Information 
    Theory, Akademiai Kiado, Budapest, pages 267-281.
R. Sleeman and T. van Eck (1999). Robust automatic P-phase picking:
    an on-line implementation in the analysis of broadband seismogram
    recordings. Phys. Earth Planet. Int., 113, 265-275, 
    doi:10.1016/S0031-9201(99)00007-2.

*************************************************************************
Author: Nils Maercklin, 6 March 2008 - 2 June 2010
*************************************************************************/
{
    register int i,k;
    float sigmanoise=0.0;
    float sigmasignal=0.0;
    float aicmax=-FLT_MAX;
    float aicmin=FLT_MAX;
    int iaicmin=0;
    int nk = n - istart;

    /* Check window length */
    if (nk<=0) return 0;

    /* Compute variances and AIC */
    for (k=1;k<nk-1;k++) {
        sigmanoise  = 0.0;
        sigmasignal = 0.0;

        /* Variances (calculation can be optimized...) */
        for (i=0;i<k;i++) {
            sigmanoise += ndata[i+istart]*ndata[i+istart];
        }
        sigmanoise /= ((float)k);

        for (i=k;i<nk-1;i++) {
            sigmasignal += sdata[i+istart]*sdata[i+istart];
        }
        sigmasignal /= ((float)(nk-k));


        /* AIC */
        if (sigmanoise>0.0 && sigmasignal>0.0) {
            aic[k+istart] = ((float)(k))*log10(sigmanoise) + \
                     ((float)(nk-k))*log10(sigmasignal) + ((float)(2*nk));

            if (aic[k+istart]>aicmax) aicmax=aic[k+istart];
            if (aic[k+istart]<aicmin) {
                aicmin  = aic[k+istart];
                iaicmin = k;
                *snr = sqrt(sigmasignal/sigmanoise);
            }
        }
        else {
            aic[k+istart]=FLT_MAX;
            *snr = 0.0;
        }
    }

    /* Fill AIC array */
    if (n>2) {
        for (i=0;i<=istart;i++) aic[i] = aicmax;
        aic[n-1] = aicmax;
    }

    /* Return index of AIC minimum (0 on failure) */
    if (iaicmin>1 && iaicmin<nk-1 && (*snr)>0.0) return iaicmin+istart;
    else return 0;
}



int aic_pick_quality(float *aic, int n, float dt, float snr, float frac, \
    float wscale, char *phase, char *pickstring, float *terr)
/************************************************************************
aic_pick_quality - generate hypo71-style pick quality string from AIC
*************************************************************************
Input:
aic         array[n] of AIC function
n           number of samples aic[]
dt          time sampling rate in seconds
snr         signal-to-noise ratio
frac        fraction of AIC range defining minimum
wscale      scale factor of hard-coded time ranges for pick quality
phase       phase letter in pick description string, e.g. P or S

Output:
            return values: 1 = pick accepted, 0 = pick rejected
pickstring  4-character pick description string (hypo71 style)
terr        estimated picking error in seconds

Note: This is an experimental function; several variables are hard-coded.
*************************************************************************
Author: Nils Maercklin, 2 June 2010 - 4 July 2011
*************************************************************************/
{
    register int i;
    int itmin, itmax, ipick;
    float tw[4];
    float trange;
    float aicmin, aicmax, aicthresh;
    int weight=4;
    char otype;

    /* Time ranges associated with weights 0-3 (adopted from ISNet) */
    tw[0] = 0.05; tw[1] = 0.1; tw[2] = 0.2; tw[3] = 0.5;

    /* Check frac and wscale */
    if (frac<0.0) frac = 0.0;
    wscale = fabs(wscale);


    /* AIC threshold */
    aicmin = aicmax = aic[0];
    ipick  = 0;
    for (i=1;i<n;i++) {
        if (aic[i]<aicmin) {
            aicmin = aic[i];
            ipick  = i;
        }
        if (aic[i]>aicmin) {
            aicmax = aic[i];
        }
    }
    aicthresh = aicmin + frac * (aicmax - aicmin);

    /* Time range for weight */
    itmin = itmax = ipick;
    /*i = 0;
    do {
        itmin = i++;
    } while (aic[i]>aicthresh);

    i = n - 1;
    do {
        itmax = i--;
    } while (aic[i]>aicthresh);*/
    for (i=ipick; i>=0; i--) {
        if (aic[i]>aicthresh) {
            itmin = i;
            break;
        }
    }
    for (i=ipick; i<n; i++) {
        if (aic[i]>aicthresh) {
            itmax = i;
            break;
        }
    }


    trange = dt * ((float)(itmax-itmin));

    /* Weight */
    for (i=0;i<4;i++) {
        if (trange<=wscale*tw[i]) {
            weight = i;
            break;
        }
    }

    /* Onset type */
    otype = (snr<5.0) ? 'E' : 'I';

    /* Pick quality string and estimated error */
    if (phase[0]=='-') phase[0]='?';
    sprintf(pickstring, "%c%c %d", otype, phase[0], weight);
    *terr = trange;

    /* Return 0, if norm. AIC range is very small; else return 1 */
    if ( (aicmax-aicmin)/((float)n) > 0.001) return 1;
    else return 0;
}



#define SQR(a) ( ((a) == 0.0) ? 0.0 : (a)*(a))
void memcof(float data[], int n, int m, float *xms, float d[])
/************************************************************************
memcof - compute linear prediction coefficients (Numerical Recipes)
*************************************************************************
Input:
data        array[n] of time series
n           number of samples in array data
m           number of output linear prediction coefficients

Output:
xms         mean square discrepancy
d           arrav[m] of linear prediction coefficients

*************************************************************************
Reference:
Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P.
    1996: Numerical Recipes in C - The Art of Scientific Computing.
    Cambridge University Press, Cambridge (Chapter 13.6).

Modifications: index ranges zero-based (NM), memory allocation (NM).
*************************************************************************/
{
    int k,j,i;
    float p=0.0,*wk1,*wk2,*wkm;

    wk1=fmalloc1(n);
    wk2=fmalloc1(n);
    wkm=fmalloc1(m);
    for (j=0;j<n;j++) p += SQR(data[j]);
    *xms=p/n;
    wk1[0]=data[0];
    wk2[n-2]=data[n-1];
    for (j=1;j<n-1;j++) {
        wk1[j]=data[j];
        wk2[j-1]=data[j];
    }
    for (k=0;k<m;k++) {
        float num=0.0,denom=0.0;
        for (j=0;j<(n-k-1);j++) {
            num += wk1[j]*wk2[j];
            denom += SQR(wk1[j])+SQR(wk2[j]);
        }
        d[k]=2.0*num/denom;
        *xms *= (1.0-SQR(d[k]));
        for (i=0;i<k;i++)
            d[i]=wkm[i]-d[k]*wkm[k-1-i];
        if (k == m-1) {
            free(wkm);
            free(wk2);
            free(wk1);
            return;
        }
        for (i=0;i<=k;i++) wkm[i]=d[i];
        for (j=0;j<(n-k-2);j++) {
            wk1[j] -= wkm[k]*wk2[j];
            wk2[j]=wk2[j+1]-wkm[k]*wk1[j+1];
        }
    }
    error("never get here in %s\n", __func__);
}
#undef SQR



void pefilter(float *data, int nt, int istart, int ncorr, int nlag, \
    float *pedata)
/************************************************************************
pefilter - apply autoregressive prediction error filter (wthout gap)
*************************************************************************
Input:
data        array[nt] of time series (seismic data)
nt          number of samples in arrays data and pedata
istart      first sample in correlation window
ncorr       number of samples in correlation window
nlag        number of lags (length of filter operator)

Output:
pedata      array[nt] of prediction-error filtered time series
            (pedata may be different from or equal to data)

*************************************************************************
Author:  Nils Maercklin, 4 March 2008  ("memcof version", 2 Sep. 2008)
*************************************************************************/
{
    register int i, j, n;
    float *fop=NULL;     /* filter coefficients (AR coefficients) */
    float xms;

    /* Allocate space for filter coefficients */
    fop   = fmalloc1(nlag);

    /* Compute autoregressive filter coefficients */
    memcof(data+istart, ncorr, nlag, &xms, fop);

    /* Apply prediction-error filter */
    /* (backward computation allows for data equal to pedata) */
    for (i=nt-1; i>=0; i--) {
        pedata[i] = data[i];
        n = (i<nlag) ? i : nlag;

        /* Convolution */
        for (j=0; j<n; j++) {
            pedata[i] -= fop[j-0] * data[i-j];
        }
    }

    /* Free space */
    free(fop);
}

/* END OF FILE */
