/* Copyright (c) 2021 Nils Maercklin */

/************************************************************************
SACAPKBK - SAC Automatic single-trace P-Phase Picker (Baer-Kradolfer)
*************************************************************************

Author: Nils Maercklin,
    RISSC, University of Naples, Italy, September 2008

Version: 2011-06-12

Version history:
    2008-03-02 (NM): Initial coding for SEG-Y traces and CWP/SU libs
    2008-09-11 (NM): Version for SAC binary traces using "nmsaclib.c"
    2008-10-06 (NM): Read multiple traces from stdin
    2011-06-12 (NM): Added flexible output file names

Notes: Input files must be evenly-sampled time series files in SAC 
    binary format (assuming NVHDR=6, IFTYPE=ITIME (1), LEVEN=1).

Reference:
M. Baer and U. Kradolfer (1987). An automatic phase picker for local 
    and teleseismic events. Bull. Seism. Soc. Am., 77(4), 1437-1445.
*************************************************************************/

#include "nmsaclib.h"     /* SAC-related functions, and other stuff */


/* Self-documentation */
const char *sdoc[] = {
"                                                                       ",
" SACAPKBK - SAC Automatic single-trace P-Phase Picker (Baer-Kradolfer) ",
"                                                                       ",
" Usage: sacapkbk [parameters] -f sac_files                             ",
"        sacapkbk < in_sac_file [parameters] > out_sac_file             ",
"                                                                       ",
" Parameters and defaults:                                              ",
"   -f        list of SAC binary waveform files (or read from stdin)    ",
"   -k    t1  SAC header field to store pick time                       ",
"   -u   1.0  minimum event uptime (CF>S1) required for trigger (sec)   ",
"   -d   u/3  maximum downtime in trigger period (seconds)              ",
"   -s    10  threshold S1 of characteristic function CF for trigger    ",
"                                                                       ",
"   -s1   10  threshold S1 of CF for trigger (same as parameter -s)     ",
"   -s2 2*s1  threshold S2 of CF for update of standard deviation       ",
"   -i   2.0  initialization time at beginning of trace (seconds)       ",
"   -e   2*i  first pick evaluation time (sec from beginning of trace)  ",
"   -p     u  phase duration for amplitude determination (seconds)      ",
"                                                                       ",
"   -b1    0  Butterworth low-cut frequency in Hz  (0 = no filter)      ",
"   -b2    0  Butterworth high-cut frequency in Hz (0 = no filter)      ",
"   -bp    3  number of poles of Butterworth filters                    ",
"                                                                       ",
"   -o        flag: overwrite input files (default appends \".out\")    ",
"   -x        flag: write filtered trace to output (flag -o ignored)    ",
"   -v        flag: verbose operation                                   ",
"                                                                       ",
"   -ext out  output file extension (if -f is set)                      ",
"   -dir      output directory (default is same as input directory)     ",
"                                                                       ",
" NM, 2011-06-12                                                        ",
"                                                                       ",
NULL};


/* Prototype of function used internally */
int ppicker_bk(float *reltrc, float *trace, int npts, int nsps, \
    int tdownmax, int tupevent, float thrshl1, float thrshl2, \
    int cf_flag, int nrpreset, int p_dur, int evaltime, \
    int *optime, float *opamp, int *opamptime, int *oifrst, \
    float *onoise, int *onoisetime, float *osignal, int *osignaltime);


int 
main(int argc, char **argv)
{
    /* Variables */
    register int i,it,iarg; /* loop indices */
    SACHEAD hd;          /* SAC header */
    float *data=NULL;    /* SAC trace data */
    float *fdata=NULL;   /* filtered SAC trace */
    float *odata=NULL;   /* pointer to output data array (e.g. data) */
    char *wkey=NULL;     /* header key, output time */
    char *outfname=NULL; /* output file name */
    const char *pol="U D";  /* polarity label */
    int isfile;          /* file names flag */
    int overwrite;       /* overwrite flag */
    int swap;            /* byte-swapping flag */
    int verbose;         /* verbose flag */
    int filtout;         /* filtered-trace output flag */
    int ntr=0;           /* trace counter (for user info) */
    int npick=0;         /* pick counter (for user info) */
    int iwkey=0;         /* header field index of wkey */

    float fmin, fmax;    /* Butterworth filter frequencies */
    float f;             /* normalized frequency */
    int npoles;          /* number of poles of Butterworth filters */

    char *outdir=NULL;   /* optional output directory name */
    char *outext=NULL;   /* output file extension */

    float tpick;     /* time of final p-phase pick in seconds */
    float dt=0.0;    /* time sampling interval in seconds */
    float ft=0.0;    /* time of first sample (delay recording time) */
    float tdmax=0.0; /* maximum downtime in trigger period (seconds) */
    float tup=0.0;   /* minimum event uptime in seconds */
    float tinit=0.0; /* initialization window length (seconds) */
    float teval=0.0; /* time of first pick evaluation (seconds) */
    float tpdur=0.0; /* phase duration for amplitude determ. (seconds) */
    float s1=0.0;    /* threshold S1 for pick flag (e.g. 10) */
    float s2=0.0;    /* threshold S2 for sigma update (e.g. 20) */
    float snr=0.0;   /* ratio of max. phase amplitude to max. noise amp. */

    float pamp=0.0;  /* maximum phase amplitude */
    float signal=0.0;/* maximum signal amplitude */
    float noise=0.0; /* maximum noise amplitude */
    int ifrst=0;     /* direction of first motion (1 or -1; 0=undetermined) */
    int ptime=0;     /* phase pick time (sample index) */
    int pamptime=0;  /* time of phase amplitude determination (sample) */
    int signaltime=0;/* time of signal amplitude (sample index) */
    int noisetime=0; /* time of noise amplitude (sample index) */
    int nsps=0;      /* number of samples per second */
    int tupevent=0;  /* minimum event uptime (number of samples) */
    int tdownmax=0;  /* maximum event downtime in trigger period (samples) */
    int nrpreset=0;  /* initialization window length (samples; BK: 2*nsps) */
    int evaltime=0;  /* time of first pick evaluation (samples; BK: 256) */
    int p_dur=0;     /* phase dur. for amp. determ. (samples; BK: 6*nsps) */


    /* Print documentation */
    if (argc==1) {
        i=0;
        while (sdoc[i]) fprintf(stderr, "%s\n", sdoc[i++]);
        return (EXIT_FAILURE);
    }

    /* Command line parameters */
    if (!fgetpar(argc, argv, "-s",  &s1))      s1=10.0;
    (void)fgetpar(argc,argv, "-s1", &s1);
    if (!fgetpar(argc, argv, "-s2", &s2))      s2=2.0*s1;
    if (!fgetpar(argc, argv, "-u",  &tup))     tup=1.0;
    if (!fgetpar(argc, argv, "-d",  &tdmax))   tdmax=tup/3.0;
    if (!fgetpar(argc, argv, "-i",  &tinit))   tinit=2.0;
    if (!fgetpar(argc, argv, "-e",  &teval))   teval=2.0*tinit;
    if (!fgetpar(argc, argv, "-p",  &tpdur))   tpdur=tup;
    if (!sgetpar(argc, argv, "-k",  &wkey))    wkey="t1";
    (void)sgetpar(argc,argv, "-wk", &wkey);

    if (!fgetpar(argc, argv, "-b1", &fmin))    fmin=0.0;
    if (!fgetpar(argc, argv, "-b2", &fmax))    fmax=0.0;
    if (!igetpar(argc, argv, "-bp", &npoles))  npoles=3;

    if (!(isfile=getflag(argc, argv, "-f")))   isfile=0;
    if (!(overwrite=getflag(argc, argv, "-o"))) overwrite=0;
    if (!(filtout=getflag(argc, argv, "-x")))  filtout=0;
    if (!(verbose=getflag(argc, argv, "-v")))  verbose=0;

    if (!sgetpar(argc, argv, "-ext", &outext)) outext="out";
    if (!sgetpar(argc, argv, "-dir", &outdir)) outdir=NULL;


    /* Get and check header field index */
    if ((iwkey=sac_hdr_index(wkey))==-1) \
        error("%s: invalid header field k=%s\n", argv[0], wkey);
    if (!IS_SACHDR_TMARK(iwkey)) \
         error("%s: header field k=%s not allowed\n", argv[0], wkey);

    /* Check other parameters */
    if (fmin<=0.0 && fmax<=0.0) npoles = 0;
    if (!isfile && isatty(1)) \
        error("%s: can't write a SAC binary file to tty\n", argv[0]);

    /* Protect original data in case of windowed output */
    if (overwrite && filtout) overwrite=0;


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

        /* Get required header values */
        dt = hd.delta;
        ft = hd.b;

        /* Check header values */
        if (ft==SAC_HEADER_FLOAT_UNDEFINED  || dt<=0.0) {
            error("%s: invalid header values: delta=%g b=%g\n", \
                argv[0], hd.delta, hd.b);
        }

        /* Convert times from seconds to samples */
        nsps     = NINT(1.0/dt);
        tupevent = NINT(tup/dt);
        tdownmax = NINT(tdmax/dt);
        nrpreset = NINT(tinit/dt);
        p_dur    = NINT(tpdur/dt);
        evaltime = NINT(teval/dt);

        /* Allocate space for additional, filtered trace */
        if (!(fdata=fmalloc1(hd.npts))) {
            error("%s: can't allocate space for trace, npts=%d\n", \
                argv[0], hd.npts);
        }
 
        /* Copy data and subtract mean */
        if (hd.depmen!=SAC_HEADER_FLOAT_UNDEFINED) {
            for (it=0; it<hd.npts; it++) fdata[it] = data[it] - hd.depmen;
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


        /* Pick P arrival time */
        ptime = ppicker_bk(fdata, data, hd.npts, nsps, tdownmax, tupevent, \
            s1, s2, 0, nrpreset, p_dur, evaltime, \
            &ptime, &pamp, &pamptime, &ifrst, &noise, &noisetime, \
            &signal, &signaltime);
        tpick = ft + dt*((float)ptime);

        /* Signal-to-noise ratio */
        snr = (ptime && noise) ? pamp/noise : 0.0;


        /* Diagnostic print */
        if (verbose) {
            fprintf(stderr, "%s: ", argv[0]);
            if (ptime) \
                fprintf(stderr, "pick %s=%g snr=%g %c\n", \
                    wkey, tpick, snr, pol[ifrst+1]);
            else fprintf(stderr, "no pick found\n");
        }

        /* Write pick time to header, increment pick counter */
        if (ptime) {
            put_sac_hdr_float(&hd, iwkey, tpick);
            npick++;
        }


        /* Select data for output (original or filtered) */
        if (filtout) {
            odata = fdata;
            set_sac_depminmax(&hd, odata);
        }
        else {
            odata = data;
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


        /* Free space for trace data */
        free(fdata);
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

int ppicker_bk(float *reltrc, float *trace, int npts, int nsps, \
    int tdownmax, int tupevent, float thrshl1, float thrshl2, \
    int cf_flag, int nrpreset, int p_dur, int evaltime, \
    int *optime, float *opamp, int *opamptime, int *oifrst, \
    float *onoise, int *onoisetime, float *osignal, int *osignaltime)
/************************************************************************
ppicker_bk - automatic P-phase picker for local and teleseismic events
             (algorithm after Baer and Kradolfer (1987), BSSA)
*************************************************************************
Input:
reltrc      array[npts] of time series, possibly filtered
trace       array[npts] of time series, used to determine amplitudes
npts        number of samples in reltrc[] and trace[]
nsps        number of samples per second
tdownmax    allowed maximum number of samples with CF<thrshl1 for pick
tupevent    required minimum number of samples with CF>thrshl1 for pick
thrshl1     threshold S1 for pick flag    (e.g. 10; BSSA paper, p. 1440)
thrshl2     threshold S2 for sigma update (e.g. 20; BSSA paper, p. 1440)
cf_flag     replace reltrc[] by characteristic function CF before return
nrpreset    time window for parameter initialization (number of samples)
p_dur       phase duration for amplitude determination (number of samples)
evaltime    time of first pick evaluation (number of samples; orig. 256)

Output:
ptime       returns pick time (sample index; ptime=0 means no pick found)
optime      pick time (sample index; same as ptime)
opamp       phase amplitude
opamptime   time of phase amplitude (sample index)
oifrst      polarity (1 or -1; 0 if no pick found)
onoise      maximum noise amplitude
onoisetime  time of maximum noise amplitude (sample index)
osignal     maximum signal amplitude
osignaltime time of maximum signal amplitude (sample index)

*************************************************************************
Notes:
This function is essentially a translation of the Fortran subroutines
in "ppicker.f" by M. Baer, keeping most of the original variable names.
Possible values for tdownmax and tupevent are given in the BSSA paper, 
pages 1440-1441. All times are passed as the sample index of the input 
trace (integer), and amplitudes are floats. 
If nrpreset, evaltime, or p_dur equals 0 in the function call, the 
values hard-coded in "ppicker.f" are used.

References:
M. Baer and U. Kradolfer (1987). An automatic phase picker for local and
    teleseismic events. Bulletin of the Seismological Society of America
    (BSSA), 77(4), 1437-1445.

The Fortran implementation by M. Baer ("ppicker.f") is available via the
ORFEUS software site: http://www.orfeus-eu.org/Software/softwarelib.html.

*************************************************************************
Author: Nils Maercklin, 2 March 2008
Credit: M. Baer, Schweizer Erdbebendienst (subroutines in "ppicker.f")
*************************************************************************/
{
    /* Variables */
    register int i;   /* sample index in loop over samples */
    int itar=0;       /* first triggered sample */
    int dtime=0;      /* counts duration where CF drops below thrshl1 */
    int amptime=0;    /* time of maximum amplitude */
    int ptime=0;      /* time of triggered phase amplitude */
    int pamptime=0;   /* time of maximum amplitude within trigger period */
    int preptime=0;   /* time of a possible earlier pick */
    int prepamptime=0;/* time of amplitude of possible earlier pick */
    int noisetime=0;  /* time of maximum noise amplitude */
    int signaltime=0; /* time of maximum signal amplitude */
    int end_dur=0;    /* end time phase amplitude determination */
    int ipkflg=0;     /* pick flag */
    int ifrst=0;      /* direction of first motion (1 or -1) */
    int uptime=0;     /* counter of samples where ipkflg>0 but no pick */
    int process=0;    /* internal processing flag */
    float ssx=0.0;    /* sum of characteristic function (CF) */
    float ssx2=0.0;   /* sum of squares of CF */
    float mean=0.0;   /* mean of CF */
    float sdev=0.0;   /* standard deviation (sigma) of CF */
    float edat=0.0;   /* function value of Equation 7 in BSSA paper */
    float edev=0.0;   /* CF value, defined in Equation 8 in BSSA paper */
    float num=0.0;    /* number of data points used for CF and mean */
    float omega=0.0;  /* weighting factor for derivative of reltrc */
    float rawold=0.0; /* last datapoint of trace when leaving function */
    float y2=0.0;     /* last value of amplitude squared */
    float yt=0.0;     /* last value of derivative squared */
    float amp=0.0;    /* maximum amplitude value */
    float ampi=0.0;   /* current amplitude (iamp in "ppicker.f") */
    float pamp=0.0;   /* maximum amplitude within trigger period */
    float prepamp=0.0;/* maximum amplitude of possible earlier pick */
    float noise=0.0;  /* maximum noise amplitude */
    float signal=0.0; /* maximum signal amplitude */
    float ysv, yy2, yyt; /* variables from subroutine "preset()" */

    int  itrm, picklength;
    float sum, rdif, rdat, rda2, rdi2;


    /* Unless positive, set parameters that are hard-coded in "ppicker.f" */
    if (nrpreset<=0) nrpreset = 2 * nsps;
    if (p_dur<=0)    p_dur    = 6 * nsps;
    if (evaltime<=0) evaltime = 256;

    /* Validate parameters */
    if (npts<=0 || nsps<=0 || tdownmax<=0 || tupevent<=tdownmax || \
        thrshl1<=0.0 || thrshl2<=0.0) {
        fprintf(stderr,"%s: invalid parameter(s)\n", __func__);
        return 0;
    }


    /* Initialize variables; subroutine "preset()" in "ppicker.f" */
    ysv    = reltrc[0];
    rawold = ysv;
    ssx    = ysv;
    y2     = 0.0;
    yt     = 0.0;
    for (i=1;i<nrpreset;i++) {
        yy2 =  reltrc[i];
        yyt =  (yy2 - ysv) * ((float)nsps);
        ysv =  yy2;
        ssx += ysv;
        y2  += yy2*yy2;
        yt  += yyt*yyt;
    }
    sdev = sqrt(((float)nrpreset)*y2-ssx*ssx)/((float)(nrpreset*nrpreset));
    ssx    = 0.0;
    ssx2   = 0.0;
    num    = 0.0;
    itar   = 0;
    ipkflg = 0;
    ptime  = 0;
    preptime = 0;

    /* Initialize remaining variables; subroutine "ppick()" in "ppicker.f" */
    ifrst  = 0;
    end_dur= 0;                              /* p_dur now in function call */
    omega  = (yt) ? y2/yt : 1.0;
    amp    = 0.0;
    pamp   = 0.0;
    prepamp= 0.0;
    noise  = 0.0;
    signal = 0.0;
    pamptime    = 0;
    prepamptime = 0;
    noisetime   = 0;
    signaltime  = 0;
    picklength  = npts - nsps;


    /* Loop over samples */
    for (i=0;i<npts;i++) {

        process = 1;                 /* flag to avoid the goto (see below) */

        if (i>picklength) {
            ampi = fabs(trace[i] + 0.5);
            if (ampi>amp) {
                amp = ampi;
                amptime = i;
            }
            process = 0;        /* "goto 160" (loop start) in "ppicker.f") */
        }
        
        if (process) {
            rdat   = reltrc[i];
            rdif   = (rdat - rawold) * ((float)nsps);
            rawold = rdat;
            rda2   = rdat*rdat;
            rdi2   = rdif*rdif;
            y2     = y2 + rda2;
            yt     = yt + rdi2;

            edat   = rda2 + omega*rdi2;        /* Equation 7 in BSSA paper */

            edat = edat*edat;
            omega = y2/yt;
            edev = (edat - mean) / sdev;       /* Equation 8 in BSSA paper */

            /* Replace reltrc by edev values (not in "ppicker.f") */
            if (cf_flag) reltrc[i] = edev;
            
            ampi = fabs(trace[i] + 0.5);

            if (ampi>amp) {
                amp = ampi;
                amptime = i;
            }
            if (i<=end_dur) {
                pamp = amp;
                pamptime = amptime;
            }

            if (edev>thrshl1 && i>evaltime) {

                if (ipkflg==0) {           /* save the current parameters  */
                    itar   = i;            /* itar is 1st triggered sample */
                    ipkflg = 1;

                    if (ptime==0) {
                        end_dur = itar + p_dur;
                        if (noise==0.0) {
                            noise = amp;
                            noisetime = amptime;
                        }
                        if (rdif<0.0) ifrst =  1;
                        if (rdif>0.0) ifrst = -1;
                    }

                    if (preptime==0) {
                        preptime = itar;
                        prepamp = amp;
                        prepamptime = amptime;
                    }
                    uptime = 1;

                }
                else if (ptime==0) {
                    if (edev>40.0 && dtime==0) ipkflg += 2;
                    uptime++;
                }

                dtime = 0;
            }
            else {
                if (ipkflg!=0) {
                    dtime++;

                    if (ptime==0) {
                        uptime++;
                    }
                    if (dtime>tdownmax) {
                        itrm = i - itar - dtime + ipkflg;

                        if (itrm>tupevent) {         /* triggered for more */
                            if (ptime==0) {          /* than tupevent      */
                                ptime = itar;
                                itar  = 0;
                            }
                        }
                        else {
                            prepamp = amp;
                            prepamptime = amptime;
                            itar = 0;
                        }

                        ipkflg = 0;
                        uptime = 0;
                    }
                }
            }

            /* Update standard deviation */
            if (edev<thrshl2 || i<=evaltime) {
                ssx  += edat;
                ssx2 += edat*edat;
                sum  =  num + 1.0;
                sdev =  sqrt( (sum*ssx2 - ssx*ssx) / (sum*sum) );

                if (sdev<=0.0) sdev = 1.0;
                mean = ssx / sum;
                num  = sum + 0.5;
            }
        } /* end process */
    } /* end loop over samples */


    /* Finalize ("if (i.gt.npts) ..." at loop start in "ppicker.f") */
    signal     = amp;
    signaltime = amptime;

    if (ptime==0 && itar!=0) {
        itrm  = i - itar - dtime + ipkflg;

        if (itrm>tupevent) {    /* triggered for more than tupevent */
            if (ptime==0) {
                ptime = itar;
                itar  = 0;
            }
        }
    }


    /* Pass variables to calling function */
    *optime      = ptime;
    *opamp       = (ptime) ? pamp     : 0.0;
    *opamptime   = (ptime) ? pamptime : 0;
    *oifrst      = (ptime) ? ifrst    : 0;
    *onoise      = noise;
    *onoisetime  = noisetime;
    *osignal     = signal;
    *osignaltime = signaltime;


    /* Return ptime (ptime=0 means no pick found) */
    return ptime;
}

/* END OF FILE */
